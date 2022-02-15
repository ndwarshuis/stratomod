import re
from os.path import basename, splitext
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import printerr
from common.bed import sort_bed_numerically

STARTCOL = "genoStart"
ENDCOL = "genoEnd"
CLASSCOL = "repClass"
FAMCOL = "repFamily"

COLS = {
    5: "genoName",
    6: STARTCOL,
    7: ENDCOL,
    11: CLASSCOL,
    12: FAMCOL,
}

RMSK_COLS = ["genoName", STARTCOL, ENDCOL]
RMSK_DF_COLS = RMSK_COLS + [CLASSCOL, FAMCOL]


def read_rmsk_df(path):
    df = read_tsv(path, header=None)[COLS].rename(columns=COLS)
    return sort_bed_numerically(df)


def merge_and_write_group(df, path, groupcol, groupname):
    dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
    merged = bt.from_dataframe(dropped).merge().to_dataframe(names=RMSK_COLS)
    if len(merged.index) == 0:
        print("WARNING: empty dataframe for %s" % path)
    else:
        merged[groupname] = merged[ENDCOL] - merged[STARTCOL]
        write_tsv(path, merged, header=True)


def parse_output(path, df, prefix):
    res = re.match("%s_(.*)" % prefix, splitext(basename(path))[0])
    if res is None:
        printerr("Unable to determine class/family from path: %s" % path)
    else:
        s = res[1].split("_")
        if len(s) == 1:
            cls = s[0]
            printerr("Filtering and merging repeat masker class %s" % cls)
            merge_and_write_group(df, path, CLASSCOL, cls)
        elif len(s) == 2:
            cls, fam = s
            printerr(
                "Filtering and merging repeat masker family %s for class %s"
                % (fam, cls)
            )
            merge_and_write_group(df, path, FAMCOL, fam)
        else:
            printerr("Invalid family/class spec in path: %s" % path)


def main():
    rmsk_df = read_rmsk_df(snakemake.input[0])
    for path in snakemake.output:
        parse_output(path, rmsk_df, snakemake.params.prefix)


main()
