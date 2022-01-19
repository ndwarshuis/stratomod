from pybedtools import BedTool as bt
from pybedtools import cleanup
from common.tsv import write_tsv
from common.cli import printerr, make_io_parser
from common.bed import sort_bed_numerically

BASE_COL = "base"
BED_COLS = ["chr", "start", "end"]
SIMPLE_REPEAT_BED_COLS = BED_COLS + [BASE_COL]

# ASSUME the input bed file for this script is sorted


def make_parser():
    parser = make_io_parser(
        "create imperfect homopolymer annotations dataframes",
        "the simple repeat bed file",
        "the merged output",
    )
    parser.add_argument(
        "-g",
        "--genome",
        required=True,
        help="the path the genome file (required to add slop)",
    )
    parser.add_argument(
        "-b",
        "--bases",
        required=True,
        help="the bases that will be merged",
    )
    return parser


def read_input(path):
    printerr("Reading dataframe from %s\n" % path)
    return bt(path).to_dataframe(names=SIMPLE_REPEAT_BED_COLS)
    # return df[df["chr"].str.match("^chr([0-9]{1,2}|X|Y)$")]


def filter_base(df, base):
    printerr("Filtering bed file for %ss" % base)
    _df = df[df[BASE_COL] == "unit=%s" % base].drop(columns=[BASE_COL])
    assert len(_df) > 0, "Filtered bed file for %s has no rows" % base
    printerr("Merging filtered bed with %i lines for %ss\n" % (len(_df), base))
    bed = bt.from_dataframe(_df).merge(d=1)
    bed.delete_temporary_history(ask=False)
    return bed


def filter_bases(df, bases):
    return [filter_base(df, b) for b in bases]


def intersect_bases(beds, bases, genome):
    slop = 5
    printerr("Intersecting merged beds for %s and adding %sbp slop" % (bases, slop))
    intersected = bt().multi_intersect(i=[b.fn for b in beds]).slop(b=slop, g=genome)
    idf = intersected.to_dataframe(header=None).iloc[:, :3]
    # these files are huge, now that we have a dataframe, remove all the bed
    # file from tmpfs to prevent a run on downloadmoreram.com
    cleanup()

    printerr("Sorting bed file")
    sortd = sort_bed_numerically(idf)

    printerr("Merging bed file")
    merged = bt.from_dataframe(sortd).merge().to_dataframe(header=None).iloc[:, :3]

    printerr("Adding homopolymer length")
    merged["%s_homopolymer_length" % bases] = merged["end"] - merged["start"] - slop * 2
    return merged


def main():
    args = make_parser().parse_args()
    # ASSUME this file is already sorted
    simreps = read_input(args.input)
    base_beds = filter_bases(simreps, args.bases)
    merged_bases = intersect_bases(base_beds, args.bases, args.genome)
    write_tsv(args.output, merged_bases, header=True)


main()
