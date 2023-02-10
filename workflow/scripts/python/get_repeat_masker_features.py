import re
import pandas as pd
from typing import Dict, Optional
import common.config as cfg
from functools import partial
from os.path import basename
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import write_tsv
from common.cli import setup_logging
from common.bed import read_bed_df

# The repeat masker database is documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

logger = setup_logging(snakemake.log[0])

# both of these columns are temporary and used to make processing easier
CLASSCOL = "_repClass"
FAMCOL = "_repFamily"

COLS = {
    11: CLASSCOL,
    12: FAMCOL,
}

fmt_feature = partial(cfg.fmt_repeat_masker_feature, snakemake.config)


def read_rmsk_df(path: str, bed_cols: Dict[str, str]):
    bed_mapping = cfg.bed_cols_indexed([5, 6, 7], bed_cols)
    prefix = cfg.refsetkey_to_chr_prefix(
        snakemake.config,
        ["annotations", "repeat_masker"],
        snakemake.wildcards["refset_key"],
    )
    return read_bed_df(path, bed_mapping, COLS, prefix, snakemake.params["filt"])


def merge_and_write_group(
    df: pd.DataFrame,
    path: str,
    bed_cols: Dict[str, str],
    groupcol: str,
    clsname: str,
    famname: Optional[str] = None,
) -> None:
    groupname = clsname if famname is None else famname
    dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
    merged = (
        bt.from_dataframe(dropped)
        .merge()
        .to_dataframe(names=cfg.bed_cols_ordered(bed_cols))
    )
    if len(merged.index) == 0:
        logger.warning("Empty dataframe for %s", path)
    else:
        col = fmt_feature(clsname, famname)
        merged[col] = merged[bed_cols["end"]] - merged[bed_cols["start"]]
        write_tsv(path, merged, header=True)


def parse_output(
    path: str,
    df: pd.DataFrame,
    file_prefix: str,
    bed_cols: Dict[str, str],
) -> None:
    res = re.match(f"{file_prefix}_(.*).tsv.gz", basename(path))
    if res is None:
        logger.error("Unable to determine class/family from path: %s", path)
    else:
        s = res[1].split("_")
        f = partial(merge_and_write_group, df, path, bed_cols)
        if len(s) == 1:
            cls = s[0]
            logger.info("Filtering and merging repeat masker class %s", cls)
            f(CLASSCOL, cls)
        elif len(s) == 2:
            cls, fam = s
            logger.info(
                "Filtering and merging repeat masker family %s for class %s",
                fam,
                cls,
            )
            f(FAMCOL, cls, fam)
        else:
            logger.info("Invalid family/class spec in path: %s", path)


def main() -> None:
    bed_cols = cfg.lookup_bed_cols(snakemake.config)
    rmsk_df = read_rmsk_df(snakemake.input[0], bed_cols)
    for path in snakemake.output:
        parse_output(path, rmsk_df, snakemake.params.file_prefix, bed_cols)


main()
