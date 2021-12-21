import sys
from functools import reduce
import pandas as pd
from pybedtools import BedTool as bt

STARTCOL = "genoStart"
ENDCOL = "genoEND"
CLASSCOL = "repClass"
RMSK_COLS = ["genoName", STARTCOL, ENDCOL]
RMSK_DF_COLS = RMSK_COLS + [CLASSCOL]


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def printerr(s):
    print(s, file=sys.stderr)


def merge_class(df, classname):
    dropped = df[df[CLASSCOL] == classname].drop(columns=[CLASSCOL])
    merged = bt.from_dataframe(dropped).merge().to_dataframe(names=RMSK_COLS)
    merged[CLASSCOL] = merged[ENDCOL] - merged[STARTCOL]
    return merged


def intersect_class(target_df, rmsk_df, cls):
    printerr("Adding annotations for class '{}'".format(cls))
    merged_bed = bt.from_dataframe(merge_class(rmsk_df, cls))
    target_cols = target_df.columns.tolist()
    return (
        bt.from_dataframe(target_df)
        .intersect(merged_bed, loj=True)
        .to_dataframe(names=target_cols + RMSK_COLS + [cls])
        .drop(columns=RMSK_COLS)
    )


def interset_rmsk(ifile, ofile, rmsk_df):
    target_df = read_tsv(ifile)
    new_df = reduce(
        lambda target, cls: intersect_class(target, rmsk_df, cls),
        ["SINE", "LINE", "LTR", "Satellite"],
        target_df,
    )
    new_df.to_csv(ofile, sep="\t", index=False, header=False)


def main():
    if len(sys.argv) != 2:
        printerr("Error: Must supply path to repeat masker tsv file")
        exit(1)

    rmsk_path = sys.argv[1]
    rmsk_df = read_tsv(rmsk_path)
    interset_rmsk(sys.stdin, sys.stdout, rmsk_df)


main()
