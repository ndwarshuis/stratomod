import sys
from functools import reduce
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser

STARTCOL = "genoStart"
ENDCOL = "genoEnd"
CLASSCOL = "repClass"
RMSK_COLS = ["genoName", STARTCOL, ENDCOL]
RMSK_DF_COLS = RMSK_COLS + [CLASSCOL]


def make_parser():
    parser = make_io_parser(
        "add repeat masker annotations to dataframe",
        "the input dataframe",
        "the annotated dataframe",
    )
    parser.add_argument(
        "-r",
        "--repeat_masker",
        required=True,
        help="the path the repeat masker TSV to use for annotations",
    )
    return parser


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
    write_tsv(ofile, new_df)


def main():
    args = make_parser().parse_args()
    interset_rmsk(args.input, args.output, read_tsv(args.repeat_masker))


main()
