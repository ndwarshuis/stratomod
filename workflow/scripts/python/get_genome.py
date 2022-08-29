from functools import partial
from common.tsv import write_tsv, read_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging
from common.functional import compose

logger = setup_logging(snakemake.log[0])

CHR_COL = 0


# NOTE only the first two columns matter here (chrom and length)
def read_input(path):
    return read_tsv(path, header=None)[[0, 1]]


def sort_genome(df):
    return df.sort_values(by=[CHR_COL], axis=0, ignore_index=True)


def write_output(path, df):
    write_tsv(path, df, header=None)


def main():
    df = compose(
        sort_genome,
        partial(standardize_chr_column, CHR_COL),
        read_input,
    )(snakemake.input[0])

    write_tsv(snakemake.output[0], df, header=None)


main()
