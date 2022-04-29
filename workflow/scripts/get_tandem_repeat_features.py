from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.bed import read_bed_df, merge_and_apply_stats, BED_START, BED_END
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


PERC_A_COL = format_base("A")
PERC_T_COL = format_base("T")
PERC_C_COL = format_base("C")
PERC_G_COL = format_base("G")

FEATURE_COLS = {
    5: "TR_unit_size",  # period
    6: "TR_unit_copies",  # copyNum
    7: "TR_consensus_size",  # consensusSize
    8: "TR_identity",  # perMatch
    9: "TR_per_indel_mismatch",  # perIndel
    10: "TR_score",  # score
    11: PERC_A_COL,  # A
    12: PERC_C_COL,  # C
    13: PERC_G_COL,  # G
    14: PERC_T_COL,  # T
}

LEN_FEATURE = "TR_length"

SLOP = 5


def read_tandem_repeats(path):
    df = read_bed_df(path, (1, 2, 3), FEATURE_COLS, snakemake.params["filt"])
    df[format_base("AT")] = df[PERC_A_COL] + df[PERC_T_COL]
    df[format_base("GC")] = df[PERC_G_COL] + df[PERC_C_COL]
    return df


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
