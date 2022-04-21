import argparse
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically

# simple wrapper script to sort bed files in a pipe using stdin/stdout

# intended usage: '... | bed_generator | \
#   THIS_SCRIPT [-H] [-c COMMENT] | \
#   bed_consumer | ...'
#
# --header: treat first line of data as the header (no header if not present)
# -c: string to use as the comment prefix
#
# NOTE: this will skip and blank/commented lines in the input


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-H", "--header", action="store_true")
    parser.add_argument("-c", "--comment", type=str, default="#")
    return parser


def main():
    args = make_parser().parse_args()
    df = read_tsv(
        None,
        header=0 if args.header is True else None,
        comment=args.comment,
        skip_blank_lines=True,
    )
    write_tsv(
        None,
        sort_bed_numerically(df, print_stderr=True),
        header=args.header,
    )


main()
