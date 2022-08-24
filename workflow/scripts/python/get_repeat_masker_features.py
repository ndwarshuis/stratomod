import re
from functools import partial
from os.path import basename, splitext
from pybedtools import BedTool as bt
from common.tsv import write_tsv
from common.cli import setup_logging
from common.bed import read_bed_df
from common.config import (
    fmt_repeat_masker_feature,
    lookup_bed_cols,
    bed_cols_indexed,
    bed_cols_ordered,
)

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

fmt_feature = partial(fmt_repeat_masker_feature, snakemake.config)


def read_rmsk_df(path, bed_cols):
    bed_mapping = bed_cols_indexed([5, 6, 7], bed_cols)
    return read_bed_df(path, bed_mapping, COLS, snakemake.params["filt"])


def merge_and_write_group(df, path, bed_cols, groupcol, clsname, famname=None):
    groupname = clsname if famname is None else famname
    dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
    merged = (
        bt.from_dataframe(dropped)
        .merge()
        .to_dataframe(names=bed_cols_ordered(bed_cols))
    )
    if len(merged.index) == 0:
        logger.warning("Empty dataframe for %s", path)
    else:
        col = fmt_feature(clsname, famname)
        merged[col] = merged[bed_cols["end"]] - merged[bed_cols["start"]]
        write_tsv(path, merged, header=True)


def parse_output(path, df, file_prefix, bed_cols):
    res = re.match(f"{file_prefix}_(.*)", splitext(basename(path))[0])
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


def main():
    bed_cols = lookup_bed_cols(snakemake.config)
    rmsk_df = read_rmsk_df(snakemake.input[0], bed_cols)
    for path in snakemake.output:
        parse_output(path, rmsk_df, snakemake.params.file_prefix, bed_cols)


main()
