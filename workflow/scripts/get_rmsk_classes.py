from pathlib import Path
import argparse
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import add_input_arg, printerr
from common.bed import sort_bed_numerically

STARTCOL = "genoStart"
ENDCOL = "genoEnd"
CLASSCOL = "repClass"
RMSK_COLS = ["genoName", STARTCOL, ENDCOL]
RMSK_DF_COLS = RMSK_COLS + [CLASSCOL]


def make_parser():
    parser = argparse.ArgumentParser(
        description="generate repeat masker TSVs for desired classes"
    )
    add_input_arg("the repeat masker source TSV", parser)
    parser.add_argument(
        "-c",
        "--classes",
        required=True,
        help="a comma-separated list of classes to filter from the source TSV",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="the directory to which to write the output",
    )
    return parser


def merge_and_write_class(df, outdir, classname):
    dropped = df[df[CLASSCOL] == classname].drop(columns=[CLASSCOL])
    merged = bt.from_dataframe(dropped).merge().to_dataframe(names=RMSK_COLS)
    merged[classname] = merged[ENDCOL] - merged[STARTCOL]
    path = Path(outdir) / ("repeat_masker_%s.tsv" % classname)
    write_tsv(path, merged, header=True)


def main():
    args = make_parser().parse_args()
    rmsk_df = sort_bed_numerically(read_tsv(args.input, names=RMSK_DF_COLS))
    for cls in args.classes.split(","):
        printerr("Filtering and merging repeat masker class %s" % cls)
        merge_and_write_class(rmsk_df, args.outdir, cls)


main()
