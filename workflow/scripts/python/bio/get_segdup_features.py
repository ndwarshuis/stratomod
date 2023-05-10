import pandas as pd
from pathlib import Path
from typing import Any, cast
import common.config as cfg
from common.tsv import write_tsv
from common.bed import read_bed, merge_and_apply_stats
from common.io import setup_logging

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
    path: Path,
    fconf: cfg.SegDupsGroup,
) -> pd.DataFrame:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = config.refsetkey_to_refkey(rsk)
    s = config.references[rk].annotations.superdups
    ocs = s.other_cols
    feature_cols = {
        ocs.align_L: fconf.fmt_col(lambda x: x.alignL),
        ocs.frac_match_indel: fconf.fmt_col(lambda x: x.fracMatchIndel),
    }
    cs = config.refsetkey_to_chr_indices(rsk)
    return read_bed(path, s.params, feature_cols, cs)


def merge_segdups(
    df: pd.DataFrame,
    fconf: cfg.SegDupsGroup,
) -> pd.DataFrame:
    bed, names = merge_and_apply_stats(fconf, df)
    return cast(pd.DataFrame, bed.to_dataframe(names=names))


def main(smk: Any, config: cfg.StratoMod) -> None:
    fconf = config.feature_names.segdups
    repeat_df = read_segdups(smk, config, smk.input[0], fconf)
    merged_df = merge_segdups(repeat_df, fconf)
    write_tsv(smk.output[0], merged_df, header=True)


# TODO make a stub so I don't need to keep repeating this
main(snakemake, snakemake.config)  # type: ignore
