import numpy as np
from typing import NamedTuple, List, Optional, Dict, Tuple, Any
import pandas as pd
import common.config as cfg
from functools import partial
from more_itertools import partition, unzip
from common.tsv import write_tsv
from common.cli import setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


class InputCol(NamedTuple):
    dtype: str
    na_value: Optional[str]


ID = "ID"
REF = "REF"
ALT = "ALT"
FORMAT = "FORMAT"
SAMPLE = "SAMPLE"


def read_vcf(input_cols: Dict[str, InputCol], path: str) -> pd.DataFrame:
    # Add columns manually because pandas can't distinguish between lines
    # starting with '##' or '#'
    dtypes, na_values = unzip(
        (
            (n, c.dtype),
            (n, c.na_value),
        )
        for n, c in input_cols.items()
    )
    df = pd.read_table(
        path,
        comment="#",
        header=None,
        # mypy complains if I don't filter on None, not sure why...
        na_values={k: v for k, v in na_values if v is not None},
        dtype={k: v for k, v in dtypes if v is not None},
        names=[*input_cols],
    )
    # NOTE: if FORMAT/SAMPLE don't exist these will just be NaNs
    return df


def fix_dot_alts(df: pd.DataFrame) -> pd.DataFrame:
    # NOTE: This step really isn't that important, but this will remove NaNs
    # from the ALT column which will make many downstream steps much easier
    df[ALT] = df[REF].where(df[ALT].isna(), df[ALT])
    return df


def lookup_maybe(k: str, d: Dict[str, float]) -> float:
    return d[k] if k in d else np.nan


def assign_format_sample_fields(
    chrom: str,
    start: str,
    fields: Dict[cfg.FeatureKey, Optional[str]],
    df: pd.DataFrame,
) -> pd.DataFrame:
    def log_split(what: str, xs: List[Tuple[cfg.FeatureKey, Optional[str]]]) -> None:
        logger.info(
            "Splitting FORMAT/SAMPLE into %i %s columns.",
            len(xs),
            what,
        )

    # Any fields that we don't want to parse will just get filled with NaNs. If
    # all fields are like this, our job is super easy and we don't need to
    # bother with zipping the FORMAT/SAMPLE columns
    fill_fields, nan_fields = map(
        lambda xs: list(xs),
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
    def split_col(n: str) -> Tuple["pd.Series[Any]", "pd.Series[int]"]:
        split_ser = df[n].str.split(":")
        len_ser = split_ser.str.len()
        return split_ser, len_ser

    format_split, format_len = split_col(FORMAT)
    sample_split, sample_len = split_col(SAMPLE)

    # warn user if any FORMAT/SAMPLE pairs are different length (which
    # shouldn't happen, but vcf files are weird)
    present = ~format_len.isna() & ~sample_len.isna()
    eqlen = format_len == sample_len

    cfg.assert_empty(
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


def get_filter_mask(
    ref_len: "pd.Series[int]",
    alt_len: "pd.Series[int]",
    filter_key: cfg.FilterKey,
) -> "pd.Series[bool]":
    snps = (ref_len == 1) & (alt_len == 1)
    keep = ~snps & ~(ref_len == alt_len) if filter_key == cfg.FilterKey.INDEL else snps
    logger.info("Number of %ss: %i", filter_key, keep.sum())
    return keep


# ASSUME there are no NaNs in alt at this point
def add_length_and_filter(
    indel_len_col: str,
    start_col: str,
    end_col: str,
    filter_key: cfg.FilterKey,
    max_ref: int,
    max_alt: int,
    df: pd.DataFrame,
) -> pd.DataFrame:
    def log_removed(
        filter_mask: "pd.Series[bool]",
        other_mask: "pd.Series[bool]",
        msg: str,
    ) -> None:
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
    len_mask = (ref_len <= max_ref) & (alt_len <= max_alt)
    log_removed(filter_mask, len_mask, f"REF > {max_ref} or ALT > {max_alt}")

    # make the output 0-based instead of 1-based (like a real bed file)
    df[start_col] = df[start_col] - 1

    df[indel_len_col] = alt_len - ref_len
    df[end_col] = df[start_col] + ref_len

    # NOTE: the main reason why all these crazy filters are in one function
    # because we get weird slice warnings unless I copy after the filter step
    # ...and for obvious reasons I don't want to copy more than once
    return df[filter_mask & len_mask & equality_mask & multi_mask].copy()


def select_columns(
    non_field_cols: List[str],
    fields: List[str],
    label_col: str,
    label: Optional[str],
    df: pd.DataFrame,
) -> pd.DataFrame:
    if label is not None:
        logger.info("Applying label %s to column %s", label, label_col)
        df[label_col] = label
    cols = non_field_cols + fields + ([] if label is None else [label_col])
    logger.info("Selecting columns for final TSV: %s", cfg.fmt_strs(cols))
    return df[cols]


def get_label(wildcards: Dict[str, str]) -> Optional[str]:
    try:
        return wildcards["label"]
    except KeyError:
        return None


# ASSUME these vcf file already have standard chromosomes
def main(smk: Any, sconf: cfg.StratoMod) -> None:
    wildcards = smk.wildcards
    fconf = sconf.feature_names
    iconf = sconf._querykey_to_input(smk.params.query_key)
    idx = fconf.bed_index

    chrom = idx.chr
    pos = idx.start
    qual = fconf.vcf.qual_name
    info = fconf.vcf.info_name
    filt = fconf.vcf.filter_name
    end = idx.end
    indel_length = fconf.vcf.len_name

    input_cols = {
        chrom: InputCol("str", None),
        pos: InputCol("int", None),
        ID: InputCol("str", "."),
        REF: InputCol("str", None),
        ALT: InputCol("str", "."),
        qual: InputCol("float", "."),
        filt: InputCol("str", "."),
        info: InputCol("str", "."),
        FORMAT: InputCol("str", "."),
        SAMPLE: InputCol("str", "."),
    }

    non_field_cols = [chrom, pos, end, indel_length, qual, filt, info]

    # TODO weakly typed
    fmt_vcf_feature = fconf.vcf.fmt_feature

    fields = {fmt_vcf_feature(k): v for k, v in iconf.format_fields.dict().items()}
    label = get_label(wildcards)

    df = read_vcf(input_cols, smk.input[0])
    # no compose operator :(
    parsed = select_columns(
        non_field_cols,
        list(fields),
        fconf.label,
        label,
        assign_format_sample_fields(
            chrom,
            pos,
            fields,
            add_length_and_filter(
                indel_length,
                pos,
                end,
                cfg.FilterKey(wildcards.filter_key),
                iconf.max_ref,
                iconf.max_alt,
                fix_dot_alts(
                    df,
                ),
            ),
        ),
    )
    return write_tsv(smk.output[0], parsed)


main(snakemake, snakemake.config)  # type: ignore
