from functools import partial
from collections import namedtuple
import numpy as np
import pandas as pd
from more_itertools import partition, unzip
from common.tsv import write_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging
from common.config import (
    fmt_vcf_feature,
    assert_empty,
    lookup_train_test_input,
    fmt_strs,
)
from common.functional import compose

logger = setup_logging(snakemake.log[0])

input_col = namedtuple("InputCol", ["dtype", "na_value"])

ID = "ID"
REF = "REF"
ALT = "ALT"
INFO = "INFO"
FORMAT = "FORMAT"
SAMPLE = "SAMPLE"


def read_vcf(input_cols, path):
    # Add columns manually because pandas can't distinguish between lines
    # starting with '##' or '#'
    dtypes, na_values = unzip(
        (
            (i, c.dtype),
            (i, c.na_value),
        )
        for i, c in enumerate(input_cols.values())
    )
    df = pd.read_table(
        path,
        comment="#",
        header=None,
        na_values=dict(na_values),
        dtype=dict(dtypes),
        names=[*input_cols],
    )
    # NOTE: if FORMAT/SAMPLE don't exist these will just be NaNs
    return df


def fix_dot_alts(df):
    # NOTE: This step really isn't that important, but this will remove NaNs
    # from the ALT column which will make many downstream steps much easier
    df[ALT] = df[REF].where(df[ALT].isna(), df[ALT])
    return df


def lookup_maybe(k, d):
    return d[k] if k in d else np.nan


def assign_format_sample_fields(chrom, start, fields, df):
    def log_split(what, xs):
        logger.info(
            "Splitting FORMAT/SAMPLE into %i %s columns.",
            len(xs),
            what,
        )

    # Any fields that we don't want to parse will just get filled with NaNs. If
    # all fields are like this, our job is super easy and we don't need to
    # bother with zipping the FORMAT/SAMPLE columns
    fill_fields, nan_fields = map(
        list,
        partition(
            lambda x: x[1] is None,
            fields.items(),
        ),
    )

    log_split("NaN", nan_fields)
    log_split("filled", fill_fields)

    for dest_col, _ in nan_fields:
        df[dest_col] = np.nan

    # If all columns are NaN'd, we are done
    if len(fill_fields) == 0:
        return df

    # For all columns that aren't to be NaN'd, split format and sample fields
    # into lists (to be zipped later).
    def split_col(n):
        split_ser = df[n].str.split(":")
        len_ser = split_ser.str.len()
        return split_ser, len_ser

    format_split, format_len = split_col(FORMAT)
    sample_split, sample_len = split_col(SAMPLE)

    # warn user if any FORMAT/SAMPLE pairs are different length (which
    # shouldn't happen, but vcf files are weird)
    present = ~format_len.isna() & ~sample_len.isna()
    eqlen = format_len == sample_len

    assert_empty(
        [
            f"{r[chrom]}@{r[start]}"
            for r in df[~eqlen & present].to_dict(orient="records")
        ],
        "Lines where FORMAT/SAMPLE have different cardinality",
    )

    # zip FORMAT/SAMPLE into dict series and assign
    #
    # NOTE: this is a faily slow operation because we need to (or at least
    # should) not assume that FORMAT/SAMPLE pair is exactly the same wrt
    # order, cardinality, etc. Some fields that are in one line might be
    # absent from others. VCF files are weird...
    mask = eqlen & present
    splits = format_split[mask].combine(
        sample_split[mask],
        lambda a, b: dict(zip(a, b)),
    )
    return pd.concat(
        [df]
        + [
            splits.apply(partial(lookup_maybe, field)).rename(dest_col)
            for dest_col, field in fill_fields
        ],
        axis=1,
    )


def get_filter_mask(ref_len, alt_len, filter_key):
    assert filter_key in ["INDEL", "SNP"], f"Invalid filter key: {filter_key}"
    snps = (ref_len == 1) & (alt_len == 1)
    keep = ~snps & ~(ref_len == alt_len) if filter_key == "INDEL" else snps
    logger.info("Number of %ss: %i", filter_key, keep.sum())
    return keep


# ASSUME there are no NaNs in alt at this point
def add_length_and_filter(
    indel_len_col,
    start_col,
    end_col,
    filter_key,
    max_ref,
    max_alt,
    df,
):
    def log_removed(filter_mask, other_mask, msg):
        logger.info(
            "Removing %i %ss with %s",
            (filter_mask & ~other_mask).sum(),
            filter_key,
            msg,
        )

    alt_len = df[ALT].str.len()
    ref_len = df[REF].str.len()
    filter_mask = get_filter_mask(ref_len, alt_len, filter_key)

    # these should be variants where ALT is a "." (which means "ALT and REF are
    # the same" and should only apply to ClinVar VCFs)
    equality_mask = df[REF] != df[ALT]
    log_removed(filter_mask, equality_mask, "same REF and ALT")

    # Remove multi-allelic variants which we cannot process (yet)
    multi_mask = ~df[ALT].str.contains(",")
    log_removed(filter_mask, multi_mask, "multi-alleles")

    # Remove any variant that doesn't pass our length filters
    len_mask = (ref_len <= max_alt) & (alt_len <= max_alt)
    log_removed(filter_mask, len_mask, f"REF > {max_ref} and LEN > {max_alt}")

    df[indel_len_col] = alt_len - ref_len
    df[end_col] = df[start_col] + ref_len
    # NOTE: the main reason why all these crazy filters are in one function
    # because we get weird slice warnings unless I copy after the filter step
    # ...and for obvious reasons I don't want to copy more than once
    return df[filter_mask & len_mask & equality_mask & multi_mask].copy()


def select_columns(non_field_cols, fields, label_col, label, df):
    if label is not None:
        logger.info("Applying label %s to column %s", label, label_col)
        df[label_col] = label
    cols = non_field_cols + [*fields] + ([] if label is None else [label_col])
    logger.info("Selecting columns for final TSV: %s", fmt_strs(cols))
    print(df)
    return df[cols]


def get_label(wildcards):
    if hasattr(wildcards, "label"):
        return wildcards.label


def main():
    wildcards = snakemake.wildcards
    sconf = snakemake.config
    fconf = sconf["features"]
    iconf = lookup_train_test_input(sconf, wildcards.input_key)
    idx = fconf["index"]

    chrom = idx["chr"]
    pos = idx["start"]
    qual = fmt_vcf_feature(sconf, "qual")
    filt = fmt_vcf_feature(sconf, "filter")
    end = idx["end"]
    indel_length = fmt_vcf_feature(sconf, "len")

    input_cols = {
        chrom: input_col(str, None),
        pos: input_col(int, None),
        ID: input_col(str, "."),
        REF: input_col(str, None),
        ALT: input_col(str, "."),
        qual: input_col(float, "."),
        filt: input_col(str, "."),
        INFO: input_col(str, "."),
        FORMAT: input_col(str, "."),
        SAMPLE: input_col(str, "."),
    }

    non_field_cols = [chrom, pos, end, indel_length, qual, filt]

    fields = {fmt_vcf_feature(sconf, k): v for k, v in iconf["format_fields"].items()}
    label = get_label(wildcards)

    return compose(
        partial(write_tsv, snakemake.output[0]),
        partial(select_columns, non_field_cols, fields, fconf["label"], label),
        partial(assign_format_sample_fields, chrom, pos, fields),
        partial(standardize_chr_column, chrom),
        partial(
            add_length_and_filter,
            fmt_vcf_feature(sconf, "len"),
            pos,
            end,
            wildcards.filter_key,
            iconf["max_ref"],
            iconf["max_alt"],
        ),
        fix_dot_alts,
        partial(read_vcf, input_cols),
    )(snakemake.input[0])


main()
