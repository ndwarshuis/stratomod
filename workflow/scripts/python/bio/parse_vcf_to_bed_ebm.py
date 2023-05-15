from typing import NamedTuple, Optional, Any, Type, Callable, cast
import pandas as pd
import common.config as cfg
from more_itertools import partition
from common.tsv import write_tsv
from common.io import setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


class InputCol(NamedTuple):
    dtype: Type[str | int | float]
    na_value: Optional[str]


ID = cfg.PandasColumn("ID")
REF = cfg.PandasColumn("REF")
ALT = cfg.PandasColumn("ALT")
FORMAT = cfg.PandasColumn("FORMAT")
SAMPLE = cfg.PandasColumn("SAMPLE")


def read_vcf(input_cols: dict[cfg.PandasColumn, InputCol], path: str) -> pd.DataFrame:
    # Add columns manually because pandas can't distinguish between lines
    # starting with '##' or '#'
    df = pd.read_table(
        path,
        comment="#",
        header=None,
        # mypy complains if I don't filter on None, not sure why...
        # na_values={k: v for k, v in na_values if v is not None},
        na_values={
            k: n for k, v in input_cols.items() if (n := v.na_value) is not None
        },
        dtype={k: v.dtype for k, v in input_cols.items() if v is not None},
        names=[*input_cols],
    )
    # NOTE: if FORMAT/SAMPLE don't exist these will just be NaNs
    return df


def fix_dot_alts(df: pd.DataFrame) -> pd.DataFrame:
    # NOTE: This step really isn't that important, but this will remove NaNs
    # from the ALT column which will make many downstream steps much easier
    df[ALT] = df[REF].where(df[ALT].isna(), df[ALT])
    return df


def assign_format_sample_fields(
    fields: dict[cfg.PandasColumn, cfg.FormatField | str | None],
    df: pd.DataFrame,
) -> pd.DataFrame:
    def log_split(what: str, n: int) -> None:
        logger.info("Splitting FORMAT/SAMPLE into %i %s columns.", n, what)

    # Any fields that we don't want to parse will just get filled with a
    # constant. If all fields are like the latter, our job is super easy and we
    # don't need to bother with zipping the FORMAT/SAMPLE columns
    fs = fields.items()
    # mypy is being dumb...
    parse_fields = [
        (k, cast(cfg.FormatField, v)) for k, v in fs if isinstance(k, cfg.FormatField)
    ]
    const_fields = [
        (k, cast(str | None, v)) for k, v in fs if not isinstance(k, cfg.FormatField)
    ]

    log_split("parsed", len(parse_fields))
    log_split("constant", len(const_fields))

    for dest_col, dest_val in const_fields:
        df[dest_col] = dest_val

    # If all columns are constant, we are done
    if len(parse_fields) == 0:
        return df

    # Else, zip FORMAT/SAMPLE into dict and assign values to column names
    #
    # NOTE: this is a fairly slow operation because we need to (or at least
    # should) not assume that FORMAT/SAMPLE pair is exactly the same wrt
    # order, cardinality, etc. Some fields that are in one line might be
    # absent from others. VCF files are weird...
    def parse(
        chrom: str,
        start: str,
        format: str,
        sample: str,
    ) -> dict[str, str | None] | tuple[str, str]:
        fs = format.split(":")
        ss = sample.split(":")
        # ASSUME any FORMAT/SAMPLE columns with different lengths are screwed
        # up in some way, so collect these to throw (semi-politely) in the
        # user's face
        if len(fs) != len(ss):
            return (chrom, start)
        d = dict(zip(fs, ss))
        return {
            col: d[field.field_name] if field.field_name in d else field.field_missing
            for col, field in parse_fields
        }

    parse_cols = [cfg.BED_CHROM, cfg.BED_START, FORMAT, SAMPLE]

    parsed = [parse(*r) for r in df[parse_cols].itertuples(index=False)]

    parsed_records, errors = partition(lambda x: isinstance(x, dict), parsed)

    if len(_errors := list(errors)) > 0:
        for chrom, start in _errors:
            logger.error(
                "FORMAT/SAMPLE have different cardinality at %s, %s", chrom, start
            )
        exit(1)

    return pd.concat([df, pd.DataFrame(parsed_records)], axis=1)


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
    df[cfg.BED_START] = df[cfg.BED_START] - 1

    df[indel_len_col] = alt_len - ref_len
    df[cfg.BED_END] = df[cfg.BED_START] + ref_len

    # NOTE: the main reason why all these crazy filters are in one function
    # because we get weird slice warnings unless I copy after the filter step
    # ...and for obvious reasons I don't want to copy more than once
    return df[filter_mask & len_mask & equality_mask & multi_mask].copy()


def select_columns(
    non_field_cols: list[cfg.PandasColumn],
    fields: list[cfg.PandasColumn],
    label_col: cfg.PandasColumn,
    label: Optional[str],
    df: pd.DataFrame,
) -> pd.DataFrame:
    if label is not None:
        logger.info("Applying label %s to column %s", label, label_col)
        df[label_col] = label
    cols = non_field_cols + fields + ([] if label is None else [label_col])
    logger.info("Selecting columns for final TSV: %s", cols)
    return df[cols]


def get_label(wildcards: dict[str, str]) -> Optional[str]:
    try:
        return wildcards["label"]
    # wildcards may look like a dict but missing keys will trigger an attr err
    except AttributeError:
        return None


# ASSUME these vcf file already have standard chromosomes
def main(smk: Any, sconf: cfg.StratoMod) -> None:
    wildcards = smk.wildcards
    fconf = sconf.feature_definitions
    iconf = sconf._querykey_to_input(smk.params.query_key)

    def fmt(f: Callable[[cfg.VCFColumns], str]) -> cfg.PandasColumn:
        return cfg.PandasColumn(fconf.vcf.fmt_name(f))

    qual = fmt(lambda x: x.qual)
    info = fmt(lambda x: x.info)
    filt = fmt(lambda x: x.filter)
    indel_length = fmt(lambda x: x.len)

    input_cols: dict[cfg.PandasColumn, InputCol] = {
        cfg.PandasColumn(cfg.BED_CHROM): InputCol(str, None),
        cfg.PandasColumn(cfg.BED_START): InputCol(int, None),
        ID: InputCol(str, "."),
        REF: InputCol(str, None),
        ALT: InputCol(str, "."),
        qual: InputCol(float, "."),
        filt: InputCol(str, "."),
        info: InputCol(str, "."),
        FORMAT: InputCol(str, "."),
        SAMPLE: InputCol(str, "."),
    }

    non_field_cols = [
        *map(cfg.PandasColumn, cfg.BED_COLS),
        indel_length,
        qual,
        filt,
        info,
    ]

    fields = {
        cfg.PandasColumn(fconf.vcf.fmt_feature(k)): v
        for k, v in iconf.format_fields.items()
    }
    label = get_label(wildcards)

    df = read_vcf(input_cols, smk.input[0])
    # no compose operator :(
    parsed = select_columns(
        non_field_cols,
        list(fields),
        fconf.label_name,
        label,
        assign_format_sample_fields(
            fields,
            add_length_and_filter(
                indel_length,
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
