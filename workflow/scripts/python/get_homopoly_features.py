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
    lookup_refset_chr_prefix,
)

logger = setup_logging(snakemake.log[0])

# temporary columns used for dataframe processing
BASE_COL = "_base"
PFCT_LEN_COL = "_perfect_length"

SLOP = 1

fmt_feature = partial(fmt_homopolymer_feature, snakemake.config)


def read_input(path, bed_cols):
    logger.info("Reading dataframe from %s", path)
    names = [*bed_cols_ordered(bed_cols), BASE_COL]
    return read_tsv(path, header=None, comment="#", names=names)


def merge_base(df, base, genome, bed_cols):
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
        .slop(b=SLOP, g=genome)
        .to_dataframe(names=[*bed_cols_ordered(bed_cols), PFCT_LEN_COL])
    )
    # these files are huge; now that we have a dataframe, remove all the bed
    # files from tmpfs to prevent a run on downloadmoreram.com
    cleanup()

    end = bed_cols["end"]
    start = bed_cols["start"]

    length_col = fmt_feature(base, "len")
    frac_col = fmt_feature(base, "imp_frac")

    merged[length_col] = merged[end] - merged[start] - SLOP * 2
    merged[frac_col] = 1 - (merged[PFCT_LEN_COL] / merged[length_col])
    return merged.drop(columns=[PFCT_LEN_COL])


def main():
    bed_cols = lookup_bed_cols(snakemake.config)
    # ASSUME this file is already sorted
    simreps = read_input(snakemake.input["bed"][0], bed_cols)
    merged = merge_base(
        simreps,
        snakemake.wildcards["base"],
        snakemake.input["genome"][0],
        bed_cols,
    )
    write_tsv(snakemake.output[0], merged, header=True)


main()
