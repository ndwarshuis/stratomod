import pandas as pd
from functools import partial
from pybedtools import BedTool as bt
from pybedtools import cleanup
from common.tsv import write_tsv, read_tsv
from common.bed import sort_bed_numerically
from common.cli import setup_logging
from common.config import (
    fmt_homopolymer_feature,
    lookup_bed_cols,
    bed_cols_ordered,
)

logger = setup_logging(snakemake.log[0])

BASE_COL = "base"

GAP_COL = "gap_count"
NORM_GAP_COL = "normalized_gap_count"
PFCT_LEN_COL = "perfect_length"

SLOP = 5

fmt_feature = partial(fmt_homopolymer_feature, snakemake.config)


def read_input(path, bed_cols):
    logger.info("Reading dataframe from %s", path)
    names = [*bed_cols_ordered(bed_cols), BASE_COL]
    return read_tsv(path, header=None, comment="#", names=names)


def filter_base(df, base, bed_cols):
    logger.info("Filtering bed file for %ss", base)
    _df = df[df[BASE_COL] == f"unit={base}"].drop(columns=[BASE_COL])
    ldf = len(_df)
    assert ldf > 0, f"Filtered bed file for {base} has no rows"
    logger.info("Merging %s rows for %ss", ldf, base)
    # Calculate the length of each "pure" homopolymer (eg just "AAAAAAAA").
    # Note that this is summed in the merge below, and the final length based
    # on start/end won't necessarily be this sum because of the -d 1 parameter
    bed_start = bed_cols["start"]
    bed_end = bed_cols["end"]
    _df[PFCT_LEN_COL] = _df[bed_end] - _df[bed_start]
    merged = (
        bt.from_dataframe(_df)
        .merge(d=1, c=[4], o=["sum"])
        .to_dataframe(names=[*bed_cols_ordered(bed_cols), PFCT_LEN_COL])
    )
    # calculate the number of "gaps" (eg imperfect homopolymer bases like the
    # "G" in "AAAAGAAAA")
    merged[GAP_COL] = merged[bed_end] - merged[bed_start] - merged[PFCT_LEN_COL]
    # these files are huge; now that we have a dataframe, remove all the bed
    # files from tmpfs to prevent a run on downloadmoreram.com
    cleanup()
    return merged


def filter_bases(df, bases, bed_cols):
    return [filter_base(df, b, bed_cols) for b in bases]


def intersect_bases(dfs, bases, genome, bed_cols):
    logger.info("Concatenating and sorting merged beds for %s", bases)
    # Assume this df is already filtered for chr1-21XY and thus we don't need
    # to do it again (save a few cpu cycles and print output)
    sdf = sort_bed_numerically(pd.concat(dfs), drop_chr=False)

    logger.info("Merging %s homopolymers and adding %sbp slop", bases, SLOP)
    # The perfect homopolymer length is column 4, and the gap length is column
    # 5; after this merge we want the sum of the perfect lengths and the max of
    # the gap lengths.
    #
    # If this is confusing, its because at some point JZ and ND had different
    # ideas of how to do this. JZ wanted the "max gap count per merged imperfect
    # homopolymer" approach and ND thought it made more sense to calculate "the
    # percentage of the final length that isn't part of a pure homopolymer."
    # (note that the latter approach considers slop and the former doesn't)
    # Regardless of which is right, they are both here ;)
    merged = (
        bt.from_dataframe(sdf)
        .slop(b=SLOP, g=genome)
        .merge(
            c=[4, 5],
            o=["sum", "max"],
        )
        .to_dataframe(
            names=[
                *bed_cols_ordered(bed_cols),
                PFCT_LEN_COL,
                GAP_COL,
            ]
        )
    )
    # put downloadmoreram.com out of business
    cleanup()

    logger.info("Adding homopolymer length/fraction features")

    length_col = fmt_feature(bases, "len")
    tot_frac_col = fmt_feature(bases, "tot_imp_frac")
    frac_col = fmt_feature(bases, "imp_frac")
    cnt_col = fmt_feature(bases, "imp_cnt")

    end = bed_cols["end"]
    start = bed_cols["start"]

    merged[length_col] = merged[end] - merged[start] - SLOP * 2
    merged[tot_frac_col] = 1 - (merged[PFCT_LEN_COL] / merged[length_col])
    merged[frac_col] = merged[GAP_COL] / merged[length_col]
    return merged.drop(columns=[PFCT_LEN_COL]).rename(columns={GAP_COL: cnt_col})


def main():
    bed_cols = lookup_bed_cols(snakemake.config)
    # ASSUME this file is already sorted
    simreps = read_input(snakemake.input["bed"][0], bed_cols)
    bases = snakemake.wildcards["bases"]
    base_dfs = filter_bases(simreps, bases, bed_cols)
    merged_bases = intersect_bases(
        base_dfs,
        bases,
        snakemake.input["genome"][0],
        bed_cols,
    )
    write_tsv(snakemake.output[0], merged_bases, header=True)


main()
