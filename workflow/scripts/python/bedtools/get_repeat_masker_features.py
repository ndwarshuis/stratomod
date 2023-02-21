import re
import pandas as pd
from typing import Optional, Any
import common.config as cfg
from os.path import basename
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import write_tsv
from common.cli import setup_logging
from common.bed import read_bed_df

# The repeat masker database is documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

logger = setup_logging(snakemake.log[0])  # type: ignore

# both of these columns are temporary and used to make processing easier
CLASSCOL = "_repClass"
FAMCOL = "_repFamily"

COLS = {
    11: CLASSCOL,
    12: FAMCOL,
}


def read_rmsk_df(smk: Any, config: cfg.StratoMod, path: str) -> pd.DataFrame:
    bed_cols = config.feature_meta.bed_index
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
    bed_cols = config.feature_meta.bed_index
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
        col = config.feature_meta.repeat_masker.fmt_name(clsname, famname)
        merged[col] = merged[bed_cols.end] - merged[bed_cols.start]
        write_tsv(path, merged, header=True)


def parse_output(config: cfg.StratoMod, path: str, df: pd.DataFrame) -> None:
    res = re.match("(.*).tsv.gz", basename(path))
    if res is None:
        logger.error("Unable to determine class/family from path: %s", path)
    else:
        s = res[1].split("_")
        if len(s) == 1:
            cls = s[0]
            logger.info("Filtering and merging repeat masker class %s", cls)
            merge_and_write_group(config, df, path, CLASSCOL, cls)
        elif len(s) == 2:
            cls, fam = s
            logger.info(
                "Filtering and merging repeat masker family %s for class %s",
                fam,
                cls,
            )
            merge_and_write_group(config, df, path, FAMCOL, cls, fam)
        else:
            logger.info("Invalid family/class spec in path: %s", path)


def main(smk: Any, config: cfg.StratoMod) -> None:
    rmsk_df = read_rmsk_df(smk, config, smk.input[0])
    for path in smk.output:
        parse_output(config, path, rmsk_df)


main(snakemake, snakemake.config)  # type: ignore
