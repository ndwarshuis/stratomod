from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.bed import read_bed_df, merge_and_apply_stats
from common.cli import setup_logging

# ASSUME segdups dataframe is fed into this script as-is. The column numbers
# below are dictionary values, and the corresponding feature names are the
# dictionary keys. Note that many feature names don't match the original column
# names in the database.
#
# CONVENTION: prepend all these features with "SEGDUP"

logger = setup_logging(snakemake.log[0])

FEATURE_COLS = {
    19: "SEGDUP_size",  # alignL
    28: "SEGDUP_identity",  # fracMatchIndel
}


def read_segdups(path):
    return read_bed_df(path, (1, 2, 3), FEATURE_COLS, snakemake.params["filt"])


def merge_segdups(df):
    stats = ["min", "max", "mean"]
    bed, names = merge_and_apply_stats(stats, "SEGDUP", df)
    return bed.to_dataframe(names=names)


def main():
    repeat_df = read_segdups(snakemake.input[0])
    merged_df = merge_segdups(repeat_df)
    write_tsv(snakemake.output[0], merged_df, header=True)


main()
