import pandas as pd
from typing import Any, cast
import common.config as cfg
from common.tsv import write_tsv
from common.bed import read_bed_df, merge_and_apply_stats
from common.cli import setup_logging

# This database is documented here:
# http://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg38&hgta_doSchemaTable=genomicSuperDups

# ASSUME segdups dataframe is fed into this script with the chromosome column
# standardized. The column numbers below are dictionary values, and the
# corresponding feature names are the dictionary keys. Note that many feature
# names don't match the original column names in the database.

logger = setup_logging(snakemake.log[0])  # type: ignore


def read_segdups(
    smk: Any,
    config: cfg.StratoMod,
    path: str,
    fconf: cfg.SegDupsGroup,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    feature_cols = {
        18: fconf.fmt_col(lambda x: x.alignL),
        27: fconf.fmt_col(lambda x: x.fracMatchIndel),
    }
    bed_mapping = bed_cols.bed_cols_indexed((1, 2, 3))
    chr_filter = config.refsetkey_to_chr_filter(
        lambda r: r.annotations.superdups.chr_prefix,
        cfg.RefsetKey(smk.wildcards["refset_key"]),
    )
    return read_bed_df(path, bed_mapping, feature_cols, chr_filter)


def merge_segdups(
    df: pd.DataFrame,
    fconf: cfg.SegDupsGroup,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    bed, names = merge_and_apply_stats(bed_cols, fconf, df)
    return cast(pd.DataFrame, bed.to_dataframe(names=names))


def main(smk: Any, config: cfg.StratoMod) -> None:
    fconf = config.feature_names.segdups
    bed_cols = config.feature_names.bed_index
    repeat_df = read_segdups(smk, config, smk.input[0], fconf, bed_cols)
    merged_df = merge_segdups(repeat_df, fconf, bed_cols)
    write_tsv(smk.output[0], merged_df, header=True)


# TODO make a stub so I don't need to keep repeating this
main(snakemake, snakemake.config)  # type: ignore
