import argparse
import logging
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically, standardize_chr_column
from common.config import refsetkey_to_chr_prefix

# simple wrapper script to sort bed files in a pipe using stdin/stdout

# intended usage: '... | bed_generator | \
#   THIS_SCRIPT [-H] [-c COMMENT] [-s CHR_COLUMN_INDEX ] | \
#   bed_consumer | ...'
#
# --header: treat first line of data as the header (no header if not present)
# -c: string to use as the comment prefix
#
# NOTE: this will skip and blank/commented lines in the input

# TODO not DRY
logging.basicConfig(level=logging.INFO)
logging.captureWarnings(True)


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-H", "--header", action="store_true")
    parser.add_argument("-c", "--comment", type=str, default="#")
    parser.add_argument("-s", "--standardize", type=int)
    parser.add_argument("-p", "--prefix", type=str, required=True)
    return parser


def main():
    args = make_parser().parse_args()
    df = read_tsv(
        None,
        header=0 if args.header is True else None,
        comment=args.comment,
        skip_blank_lines=True,
    )
    if args.standardize is not None:
        cols = df.columns.tolist()
        df = standardize_chr_column(args.prefix, cols[args.standardize], df)
    write_tsv(
        None,
        sort_bed_numerically(df),
        header=args.header,
    )


main()
