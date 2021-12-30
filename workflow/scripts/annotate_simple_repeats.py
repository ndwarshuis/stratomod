import sys
from itertools import repeat, chain
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser


MERGE_STATS = ["max", "min", "median"]
COL_INDICES = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
BED_HEADERS = ["_chr", "_chromStart", "_chromEnd"]
LEN_FEATURE = "region_length"

# change this to select different features
FEATURES = ["period_min", "copyNum_max", "perMatch_median", LEN_FEATURE]


def make_parser():
    parser = make_io_parser(
        "add simple repeat annotations to dataframe",
        "the input dataframe",
        "the annotated dataframe",
    )
    parser.add_argument(
        "-s",
        "--simple_reps",
        required=True,
        help="the path the simple repeats TSV to use for annotations",
    )
    parser.add_argument(
        "-g",
        "--genome_file",
        required=True,
        help="the path the genome file (required to add slop)",
    )
    return parser


def repeat_opt(o):
    return list(repeat(o, len(COL_INDICES)))


def printerr(s):
    print(s, file=sys.stderr)


def read_simple_repeats(path):
    printerr("Reading simple repeats")
    return read_tsv(path).drop(columns=["bin"])


def merge_simple_repeats(gfile, df):
    printerr("Merging repeat regions")

    opts = list(chain(*map(repeat_opt, MERGE_STATS)))
    cols = list(chain(*repeat(COL_INDICES, len(MERGE_STATS))))
    bed_cols = df.columns.tolist()

    # just use one column for count since all columns will produce the same
    # number
    mopts = ["count"] + opts
    mcols = [5] + cols
    headers = (
        BED_HEADERS
        + ["count"]
        + ["{}_{}".format(bed_cols[c - 1], o) for c, o in zip(cols, opts)]
    )

    merged_df = (
        bt.from_dataframe(df)
        .merge(c=mcols, o=mopts)
        .slop(b=5, g=gfile)
        .to_dataframe(names=headers)
    )
    # use the original dataframe to get the region length since we added slop to
    # the merged version
    merged_df[LEN_FEATURE] = df["chromEnd"] - df["chromStart"]

    return merged_df


def intersect_repeats(ifile, ofile, repeats_df):
    printerr("Intersecting repeats")

    target_df = read_tsv(ifile)

    intersect_headers = target_df.columns.tolist() + repeats_df.columns.tolist()

    out_df = (
        bt.from_dataframe(target_df)
        .intersect(bt.from_dataframe(repeats_df), loj=True)
        .to_dataframe(names=intersect_headers)
        .drop(columns=BED_HEADERS)
    )
    write_tsv(ofile, out_df)


def main():
    args = make_parser().parse_args()
    repeat_df = read_simple_repeats(args.simple_reps)
    merged_repeat_df = merge_simple_repeats(args.genome_file, repeat_df)
    intersect_repeats(
        args.input,
        args.output,
        merged_repeat_df[BED_HEADERS + FEATURES],
    )


main()
