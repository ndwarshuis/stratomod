from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.bed import (
    sort_bed_numerically,
    filter_chromosomes,
    merge_and_apply_stats,
    BED_CHR,
    BED_START,
    BED_END,
)
from common.cli import setup_logging

logger = setup_logging(snakemake.log[0])

# ASSUME segdups dataframe is fed into this script as-is. The column numbers
# below are dictionary values, and the corresponding feature names are the
# dictionary keys. Note that many feature names don't match the original column
# names in the database.
#
# CONVENTION: prepend all these features with "segdup"

BED_COLUMNS = {
    BED_CHR: 1,  # chrom
    BED_START: 2,  # chromStart
    BED_END: 3,  # chromEnd
}

FEATURE_COLUMNS = {
    "segdup_size": 19,  # alignL
    "segdup_identity": 28,  # fracMatchIndel
}


ALL_COLUMNS = {**BED_COLUMNS, **FEATURE_COLUMNS}


def read_segdups(path):
    df = read_tsv(path, header=None)[[*ALL_COLUMNS.values()]]
    df.columns = [*ALL_COLUMNS]
    df = filter_chromosomes(df, snakemake.params["filt"])
    return sort_bed_numerically(df)


def merge_segdups(df):
    stats = ["min", "max", "mean"]
    bed, names = merge_and_apply_stats(stats, "segdup", df)
    return bed.to_dataframe(names=names)


def main():
    repeat_df = read_segdups(snakemake.input[0])
    merged_df = merge_segdups(repeat_df)
    write_tsv(snakemake.output[0], merged_df, header=True)


main()
