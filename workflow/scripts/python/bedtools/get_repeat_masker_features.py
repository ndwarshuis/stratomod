import pandas as pd
from typing import Optional, Any
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import write_tsv
from common.io import setup_logging
from common.bed import read_bed_df

# The repeat masker database is documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

logger = setup_logging(snakemake.log[0])  # type: ignore

# both of these columns are temporary and used to make processing easier
CLASSCOL = "_repClass"
FAMCOL = "_repFamily"

COLS = {
    11: cfg.PandasColumn(CLASSCOL),
    12: cfg.PandasColumn(FAMCOL),
}


def read_rmsk_df(smk: Any, config: cfg.StratoMod, path: str) -> pd.DataFrame:
    bed_cols = config.feature_names.bed_index
    bed_mapping = bed_cols.bed_cols_indexed((5, 6, 7))
    chr_filter = config.refsetkey_to_chr_filter(
        lambda r: r.annotations.superdups.chr_prefix,
        cfg.RefsetKey(smk.wildcards["refset_key"]),
    )
    return read_bed_df(path, bed_mapping, COLS, chr_filter)


def merge_and_write_group(
    config: cfg.StratoMod,
    df: pd.DataFrame,
    path: str,
    groupcol: str,
    clsname: str,
    famname: Optional[str] = None,
) -> None:
    bed_cols = config.feature_names.bed_index
    groupname = clsname if famname is None else famname
    dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
    merged = (
        bt.from_dataframe(dropped)
        .merge()
        .to_dataframe(names=bed_cols.bed_cols_ordered())
    )
    if len(merged.index) == 0:
        logger.warning("Empty dataframe for %s", path)
    else:
        col = config.feature_names.repeat_masker.fmt_name(clsname, famname)
        merged[col] = merged[bed_cols.end] - merged[bed_cols.start]
        write_tsv(path, merged, header=True)


def parse_output(config: cfg.StratoMod, path: str, key: str, df: pd.DataFrame) -> None:
    s = key.split("_")
    if len(s) == 1:
        cls = s[0]
        logger.info("Filtering/merging rmsk class %s", cls)
        merge_and_write_group(config, df, path, CLASSCOL, cls)
    elif len(s) == 2:
        cls, fam = s
        logger.info("Filtering/merging rmsk family %s/class %s", fam, cls)
        merge_and_write_group(config, df, path, FAMCOL, cls, fam)
    else:
        logger.info("Invalid family/class spec in path: %s", path)


def main(smk: Any, config: cfg.StratoMod) -> None:
    rmsk_df = read_rmsk_df(smk, config, smk.input[0])
    for key, path in smk.output.items():
        parse_output(config, path, key, rmsk_df)


main(snakemake, snakemake.config)  # type: ignore
