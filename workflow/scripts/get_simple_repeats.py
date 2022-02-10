from itertools import repeat, chain, product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser, printerr
from common.bed import sort_bed_numerically

# columns from here: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
# ASSUME the input dataframe is the table downloaded as-is from the above src

BED_START = "chromStart"
BED_END = "chromEnd"

BED_COLUMNS = {
    "chrom": 1,
    BED_START: 2,
    BED_END: 3,
}

FEATURE_COLUMNS = {
    "period": 5,
    "copyNum": 6,
    "consensusSize": 7,
    "perMatch": 8,
    "perIndel": 9,
    "score": 10,
    "A": 11,
    "C": 12,
    "G": 13,
    "T": 14,
}

ALL_COLUMNS = {**BED_COLUMNS, **FEATURE_COLUMNS}

MERGE_STATS = ["max", "min", "median"]
LEN_FEATURE = "region_length"


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


def read_simple_repeats(path):
    df = read_tsv(path, header=None)[list(ALL_COLUMNS.values())]
    df.columns = list(ALL_COLUMNS)
    df["AT"] = df["A"] + df["T"]
    df["GC"] = df["G"] + df["C"]
    return sort_bed_numerically(df)


def merge_simple_repeats(gfile, df):
    # compute stats on all columns except the first 3
    drop_n = 3
    stat_cols = df.columns.tolist()[drop_n:]

    printerr("Computing stats for columns: %s\n" % ", ".join(stat_cols))
    printerr("Stats to compute: %s\n" % ", ".join(MERGE_STATS))

    cols, opts, headers = unzip(
        (i + drop_n + 1, m, "%s_%s" % (s, m))
        for (i, s), m in product(enumerate(stat_cols), MERGE_STATS)
    )

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ["count"] + list(opts)
    full_cols = [drop_n + 1] + list(cols)
    full_headers = list(BED_COLUMNS) + ["count"] + list(headers)

    printerr("Merging repeat regions.")
    # TODO there might be a way to make pybedtools echo what it is doing, but
    # for now this is a sanity check that this crazy command is executed
    # correctly
    printerr(
        "Using command: 'bedtools merge -i <file> -c %s -o %s'"
        % (", ".join(map(str, full_cols)), ", ".join(full_opts))
    )

    merged_df = (
        bt.from_dataframe(df)
        .merge(c=full_cols, o=full_opts)
        .slop(b=5, g=gfile)
        .to_dataframe(names=full_headers)
    )
    # use the original dataframe to get the region length since we added slop to
    # the merged version
    merged_df[LEN_FEATURE] = df[BED_END] - df[BED_START]

    return merged_df


def main():
    args = make_parser().parse_args()
    repeat_df = read_simple_repeats(args.input)
    merged_repeat_df = merge_simple_repeats(args.genome_file, repeat_df)
    write_tsv(args.output, merged_repeat_df, header=True)


main()
