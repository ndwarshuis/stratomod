import pandas as pd
from typing import Any
import common.config as cfg
from common.tsv import write_tsv
from common.bed import read_bed_df, merge_and_apply_stats
from common.cli import setup_logging

# Input dataframe documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
#
# ASSUME this dataframe is fed into this script as-is. The column numbers below
# are dictionary values, and the corresponding feature names are the dictionary
# keys. Note that many feature names don't match the original column names in
# the database.

logger = setup_logging(snakemake.log[0])  # type: ignore

SLOP = 5


def format_base(bs_prefix: int, base: str) -> str:
    return f"{bs_prefix}_{base}"


def read_tandem_repeats(
    smk: Any,
    path: str,
    fconf: cfg.TandemRepeatMeta,
    bed_cols: cfg.BedIndex,
    sconf: cfg.StratoMod,
) -> pd.DataFrame:
    fmt_base = fconf.fmt_name_base
    cols = fconf.columns
    perc_a_col = fmt_base("A")
    perc_t_col = fmt_base("T")
    perc_c_col = fmt_base("C")
    perc_g_col = fmt_base("G")
    feature_cols = {
        5: cols.period,
        6: cols.copyNum,
        8: cols.perMatch,
        9: cols.perIndel,
        10: cols.score,
        11: perc_a_col,
        12: perc_c_col,
        13: perc_g_col,
        14: perc_t_col,
    }
    bed_mapping = bed_cols.bed_cols_indexed((1, 2, 3))
    chr_filter = sconf.refsetkey_to_chr_filter(
        lambda r: r.annotations.simreps.chr_prefix,
        cfg.RefsetKey(smk.wildcards["refset_key"]),
    )
    df = read_bed_df(path, bed_mapping, feature_cols, chr_filter)
    df[fmt_base("AT")] = df[perc_a_col] + df[perc_t_col]
    df[fmt_base("GC")] = df[perc_g_col] + df[perc_c_col]
    # Filter out all TRs that have period == 1, since those by definition are
    # homopolymers. NOTE, there is a difference between period and consensusSize
    # in this database; however, it turns out that at least for GRCh38 that the
    # sets of TRs where either == 1 are identical, so just use period here
    # since I can easily refer to it.
    logger.info("Removing TRs with unitsize == 1")
    return df[df[cols.period] > 1]


def merge_tandem_repeats(
    gfile: str,
    df: pd.DataFrame,
    fconf: cfg.TandemRepeatMeta,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    prefix = fconf.prefix
    bed, names = merge_and_apply_stats(bed_cols, fconf, df)
    merged_df = bed.slop(b=SLOP, g=gfile).to_dataframe(names=names)
    len_col = f"{prefix}_{fconf.other.len}"
    merged_df[len_col] = merged_df[bed_cols.end] - merged_df[bed_cols.start] - SLOP * 2
    return merged_df


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    i = smk.input
    bed_cols = sconf.feature_meta.bed_index
    fconf = sconf.feature_meta.tandem_repeats
    repeat_df = read_tandem_repeats(smk, i.src[0], fconf, bed_cols, sconf)
    merged_df = merge_tandem_repeats(i.genome[0], repeat_df, fconf, bed_cols)
    write_tsv(smk.output[0], merged_df, header=True)


main(snakemake, snakemake.config)  # type: ignore
