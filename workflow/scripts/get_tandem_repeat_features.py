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

# Input dataframe from here: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
#
# ASSUME this dataframe is fed into this script as-is. The column numbers below
# are dictionary values, and the corresponding feature names are the dictionary
# keys. Note that many feature names don't match the original column names in
# the database.
#
# CONVENTION: prepend all these features with "TR" ("tandem repeat")


def format_base(base):
    return f"TR_percent_{base}"


BED_COLUMNS = {
    BED_CHR: 1,  # chrom
    BED_START: 2,  # chromStart
    BED_END: 3,  # chromEnd
}

PERC_A_COL = format_base("A")
PERC_T_COL = format_base("T")
PERC_C_COL = format_base("C")
PERC_G_COL = format_base("G")

FEATURE_COLUMNS = {
    "TR_unit_size": 5,  # period
    "TR_unit_copies": 6,  # copyNum
    "TR_consensus_size": 7,  # consensusSize
    "TR_identity": 8,  # perMatch
    "TR_per_indel_mismatch": 9,  # perIndel
    "TR_score": 10,  # score
    PERC_A_COL: 11,  # A
    PERC_C_COL: 12,  # C
    PERC_G_COL: 13,  # G
    PERC_T_COL: 14,  # T
}

LEN_FEATURE = "TR_length"

ALL_COLUMNS = {**BED_COLUMNS, **FEATURE_COLUMNS}

SLOP = 5


def read_tandem_repeats(path):
    df = read_tsv(path, header=None)[[*ALL_COLUMNS.values()]]
    df.columns = [*ALL_COLUMNS]
    df = filter_chromosomes(df, snakemake.params["filt"])
    df[format_base("AT")] = df[PERC_A_COL] + df[PERC_T_COL]
    df[format_base("GC")] = df[PERC_G_COL] + df[PERC_C_COL]
    return sort_bed_numerically(df)


def merge_tandem_repeats(gfile, df):
    stats = ["max", "min", "median"]
    bed, names = merge_and_apply_stats(stats, "TR", df)

    merged_df = bed.slop(b=SLOP, g=gfile).to_dataframe(names=names)
    # use the original dataframe to get the region length since we added slop to
    # the merged version
    merged_df[LEN_FEATURE] = df[BED_END] - df[BED_START]

    return merged_df


def main():
    i = snakemake.input
    repeat_df = read_tandem_repeats(i.src[0])
    merged_df = merge_tandem_repeats(i.genome[0], repeat_df)
    write_tsv(snakemake.output[0], merged_df, header=True)


main()
