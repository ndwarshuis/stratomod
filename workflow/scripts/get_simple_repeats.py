from itertools import repeat, chain
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser, printerr


MERGE_STATS = ["max", "min", "median"]
COL_INDICES = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
BED_HEADERS = ["chr", "chromStart", "chromEnd"]
LEN_FEATURE = "region_length"

# change this to select different features
FEATURES = ["period_min", "copyNum_max", "perMatch_median", LEN_FEATURE]


def make_parser():
    parser = make_io_parser(
        "merge simple repeat annotations dataframe",
        "the source dataframe",
        "the merged dataframe",
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


def read_simple_repeats(path):
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


def main():
    args = make_parser().parse_args()
    repeat_df = read_simple_repeats(args.input)
    merged_repeat_df = merge_simple_repeats(args.genome_file, repeat_df)
    write_tsv(args.output, merged_repeat_df)


main()
