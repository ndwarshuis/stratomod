import argparse
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically

# simple wrapper script to sort bed files in a pipe using stdin/stdout

# intended usage: '... | bed_generator | THIS_SCRIPT [-h] | bed_consumer | ...'
#
# -h: treat first line of data as the header (no header if not present)
#
# NOTE: this will skip and blank/commented lines in the input


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", action="store_true")
    return parser


def main():
    args = make_parser().parse_args()
    df = read_tsv(
        None,
        header=0 if args.header is True else None,
        skip_blank_lines=True,
    )
    write_tsv(None, sort_bed_numerically(df), header=args.header)


main()
