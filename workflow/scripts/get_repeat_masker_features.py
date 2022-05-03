import re
from os.path import basename, splitext
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.bed import read_bed_df, BED_CHR, BED_START, BED_END

# CONVENTION: prefix all features with "REPMASK"

logger = setup_logging(snakemake.log[0])

# both of these columns are temporary and used to make processing easier
CLASSCOL = "_repClass"
FAMCOL = "_repFamily"

COLS = {
    11: CLASSCOL,
    12: FAMCOL,
}

RMSK_COLS = [BED_CHR, BED_START, BED_END]


def read_rmsk_df(path):
    return read_bed_df(path, (5, 6, 7), COLS, snakemake.params["filt"])


def merge_and_write_group(df, path, groupcol, groupname):
    dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
    merged = bt.from_dataframe(dropped).merge().to_dataframe(names=RMSK_COLS)
    if len(merged.index) == 0:
        logger.warning("Empty dataframe for %s", path)
    else:
        merged[f"REPMASK_{groupname}_length"] = merged[BED_END] - merged[BED_START]
        write_tsv(path, merged, header=True)


def parse_output(path, df, prefix):
    res = re.match(f"{prefix}_(.*)", splitext(basename(path))[0])
    if res is None:
        logger.error("Unable to determine class/family from path: %s", path)
    else:
        s = res[1].split("_")
        if len(s) == 1:
            cls = s[0]
            logger.info("Filtering and merging repeat masker class %s", cls)
            merge_and_write_group(df, path, CLASSCOL, cls)
        elif len(s) == 2:
            cls, fam = s
            logger.info(
                "Filtering and merging repeat masker family %s for class %s",
                fam,
                cls,
            )
            merge_and_write_group(df, path, FAMCOL, fam)
        else:
            logger.info("Invalid family/class spec in path: %s", path)


def main():
    rmsk_df = read_rmsk_df(snakemake.input[0])
    for path in snakemake.output:
        parse_output(path, rmsk_df, snakemake.params.prefix)


main()
