from functools import reduce
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser, printerr


def make_parser():
    parser = make_io_parser(
        "add annotations to dataframe",
        "the input dataframe",
        "the annotated dataframe",
    )
    parser.add_argument(
        "-t",
        "--tsvs",
        required=True,
        nargs="+",
        help="the paths to the TSV files to add as annotations (space separated)",
    )
    return parser


def left_outer_intersect(left, path):
    # Use bedtools to perform left-outer join of two bed/tsv files. Since
    # bedtools will join all columns from the two input files, keep track of the
    # width of the left input file so that the first three columns of the right
    # input (chr, chrStart, chrEnd, which are redundant) can be dropped.
    left_cols = left.columns.tolist()
    left_width = len(left_cols)
    right = read_tsv(path)
    right_cols = right.columns.tolist()
    right_bed = bt.from_dataframe(right)
    printerr("Adding annotations from %s" % path)
    printerr("Annotations added: %s\n" % ", ".join(right_cols[3:]))
    new_df = (
        bt.from_dataframe(left)
        .intersect(right_bed, loj=True)
        .to_dataframe(names=left_cols + right_cols)
    )
    return new_df.drop(columns=new_df.columns[left_width : left_width + 3])


def intersect_tsvs(ifile, ofile, tsv_paths):
    target_df = read_tsv(ifile)
    new_df = reduce(left_outer_intersect, tsv_paths, target_df)
    write_tsv(ofile, new_df)


def main():
    args = make_parser().parse_args()
    printerr("Adding annotations to %s\n" % args.input)
    intersect_tsvs(args.input, args.output, args.tsvs)


main()
