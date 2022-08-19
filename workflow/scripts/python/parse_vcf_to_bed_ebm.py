from functools import partial
from collections import namedtuple
import numpy as np
import pandas as pd
from more_itertools import partition, unzip
from common.tsv import write_tsv
from common.bed import standardize_chr_series
from common.cli import setup_logging
from common.config import fmt_vcf_feature, assert_empty, lookup_train_test_input
from common.prepare import compose

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


def standardize_chrs(chr_col, df):
    df[chr_col] = standardize_chr_series(df[chr_col])
    return df.dropna(subset=[chr_col]).astype({chr_col: int})


def fix_dot_alts(df):
    # TODO not sure if these should be filled or outright removed
    df[ALT] = df[REF].where(df[ALT].isna(), df[ALT])
    return df


def remove_multialleles(df):
    # ASSUME there are no NaNs in alt at this point
    return df[~df[ALT].str.contains(",")]


def lookup_maybe(k, d):
    return d[k] if k in d else np.nan


def assign_format_sample_fields(chrom, start, fields, df):
    fill_fields_gen, nan_fields_gen = partition(
        lambda x: x[1] is None,
        fields.items(),
    )
    for dest_col, _ in nan_fields_gen:
        df[dest_col] = np.nan

    fill_fields = [*fill_fields_gen]

    # Any fields that we don't want to parse will just get filled with NaNs. If
    # all fields are like this, our job is super easy and we don't need to
    # bother with zipping the FORMAT/SAMPLE columns
    if len(fill_fields) == 0:
        return df

    # split format and sample fields into lists
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


def filter_mask(ref_len, alt_len, filter_key):
    assert filter_key in ["INDEL", "SNP"], f"Invalid filter key: {filter_key}"
    snps = (ref_len == 1) & (alt_len == 1)
    return ~snps & ~(ref_len == alt_len) if filter_key == "INDEL" else snps


def add_length_and_filter(indel_len_col, start_col, end_col, filter_key, df):
    # ASSUME there are no NaNs in alt at this point
    multi_alts = df[ALT].str.contains(",")
    alt_len = df[ALT].str.len()
    ref_len = df[REF].str.len()
    df[indel_len_col] = alt_len - alt_len
    df[end_col] = df[start_col] + ref_len
    mask = filter_mask(ref_len, alt_len, filter_key) & ~multi_alts
    return df[mask].copy()


def select_nonlabel_columns(non_field_cols, fields, df):
    return df[non_field_cols + [*fields]]


def get_label(wildcards):
    if hasattr(wildcards, "label"):
        return wildcards.label


def add_label_maybe(label, label_col, df):
    if label is not None:
        df[label_col] = label
    return df


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
        partial(add_label_maybe, label, fconf["label"]),
        partial(select_nonlabel_columns, non_field_cols, fields),
        partial(assign_format_sample_fields, chrom, pos, fields),
        partial(standardize_chrs, chrom),
        partial(
            add_length_and_filter,
            fmt_vcf_feature(sconf, "len"),
            pos,
            end,
            wildcards.filter_key,
        ),
        remove_multialleles,
        fix_dot_alts,
        partial(read_vcf, input_cols),
    )(snakemake.input[0])


main()
