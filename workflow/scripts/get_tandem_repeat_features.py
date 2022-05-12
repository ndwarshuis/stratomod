from functools import partial
from pybedtools import BedTool as bt
from common.tsv import write_tsv
from common.bed import read_bed_df, merge_and_apply_stats
from common.cli import setup_logging
from common.config import (
    fmt_tandem_repeat_base,
    lookup_bed_cols,
    bed_cols_indexed,
)

# Input dataframe documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
#
# ASSUME this dataframe is fed into this script as-is. The column numbers below
# are dictionary values, and the corresponding feature names are the dictionary
# keys. Note that many feature names don't match the original column names in
# the database.

logger = setup_logging(snakemake.log[0])

SLOP = 5


def format_base(bs_prefix, base):
    return f"{bs_prefix}_{base}"


def read_tandem_repeats(path, fconf, bed_cols):
    fmt_base = partial(fmt_tandem_repeat_base, snakemake.config)
    cols = fconf["columns"]
    perc_a_col = fmt_base("A")
    perc_t_col = fmt_base("T")
    perc_c_col = fmt_base("C")
    perc_g_col = fmt_base("G")
    feature_cols = {
        5: cols["period"],
        6: cols["copyNum"],
        8: cols["perMatch"],
        9: cols["perIndel"],
        10: cols["score"],
        11: perc_a_col,
        12: perc_c_col,
        13: perc_g_col,
        14: perc_t_col,
    }
    bed_mapping = bed_cols_indexed([1, 2, 3], bed_cols)
    df = read_bed_df(path, bed_mapping, feature_cols, snakemake.params["filt"])
    df[fmt_base("AT")] = df[perc_a_col] + df[perc_t_col]
    df[fmt_base("GC")] = df[perc_g_col] + df[perc_c_col]
    return df


def merge_tandem_repeats(gfile, df, fconf, bed_cols):
    prefix = fconf["prefix"]
    bed, names = merge_and_apply_stats(
        fconf["operations"],
        bed_cols,
        prefix,
        df,
    )

    merged_df = bed.slop(b=SLOP, g=gfile).to_dataframe(names=names)
    # use the original dataframe to get the region length since we added slop to
    # the merged version
    len_col = f"{prefix}_{fconf['other']['len']}"
    merged_df[len_col] = df[bed_cols["end"]] - df[bed_cols["start"]]

    return merged_df


def main():
    i = snakemake.input
    bed_cols = lookup_bed_cols(snakemake.config)
    fconf = snakemake.config["features"]["tandem_repeats"]
    repeat_df = read_tandem_repeats(i.src[0], fconf, bed_cols)
    merged_df = merge_tandem_repeats(i.genome[0], repeat_df, fconf, bed_cols)
    write_tsv(snakemake.output[0], merged_df, header=True)


main()
