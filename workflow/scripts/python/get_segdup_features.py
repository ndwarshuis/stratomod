import pandas as pd
from typing import Any
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
    fconf: cfg.SegDupsMeta,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    cols = fconf.columns
    feature_cols = {
        18: cols.alignL,
        27: cols.fracMatchIndel,
    }
    bed_mapping = cfg.bed_cols_indexed([1, 2, 3], bed_cols)
    prefix = cfg.refsetkey_to_ref(
        config,
        smk.wildcards["refset_key"],
    ).annotations.superdups.chr_prefix
    return read_bed_df(
        path,
        bed_mapping,
        feature_cols,
        prefix,
        smk.params["filt"],
    )


def merge_segdups(
    df: pd.DataFrame,
    fconf: cfg.SegDupsMeta,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    bed, names = merge_and_apply_stats(
        list(fconf.operations),
        bed_cols,
        fconf.prefix,
        df,
    )
    return bed.to_dataframe(names=names)


def main(smk, config: cfg.StratoMod) -> None:
    fconf = config.feature_meta.segdups
    bed_cols = cfg.lookup_bed_cols(config)
    repeat_df = read_segdups(smk, config, smk.input[0], fconf, bed_cols)
    merged_df = merge_segdups(repeat_df, fconf, bed_cols)
    write_tsv(smk.output[0], merged_df, header=True)


# TODO make a stub so I don't need to keep repeating this
main(snakemake, snakemake.config)  # type: ignore
