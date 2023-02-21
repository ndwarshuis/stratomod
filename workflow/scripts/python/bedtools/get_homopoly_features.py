import pandas as pd
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from pybedtools import cleanup
from common.tsv import write_tsv, read_tsv
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])  # type: ignore

# temporary columns used for dataframe processing
BASE_COL = "_base"
PFCT_LEN_COL = "_perfect_length"

SLOP = 1


def read_input(path: str, bed_cols: cfg.BedIndex):
    logger.info("Reading dataframe from %s", path)
    names = [*bed_cols.bed_cols_ordered(), BASE_COL]
    return read_tsv(path, header=None, comment="#", names=names)


def merge_base(
    config: cfg.StratoMod,
    df: pd.DataFrame,
    base: cfg.Base,
    genome: str,
    bed_cols: cfg.BedIndex,
) -> pd.DataFrame:
    logger.info("Filtering bed file for %ss", base)
    _df = df[df[BASE_COL] == f"unit={base.value}"].drop(columns=[BASE_COL])
    ldf = len(_df)
    assert ldf > 0, f"Filtered bed file for {base} has no rows"
    logger.info("Merging %s rows for %ss", ldf, base)
    # Calculate the length of each "pure" homopolymer (eg just "AAAAAAAA").
    # Note that this is summed in the merge below, and the final length based
    # on start/end won't necessarily be this sum because of the -d 1 parameter
    bed_start = bed_cols.start
    bed_end = bed_cols.end
    _df[PFCT_LEN_COL] = _df[bed_end] - _df[bed_start]
    merged = (
        bt.from_dataframe(_df)
        .merge(d=1, c=[4], o=["sum"])
        .slop(b=SLOP, g=genome)
        .to_dataframe(names=[*bed_cols.bed_cols_ordered(), PFCT_LEN_COL])
    )
    # these files are huge; now that we have a dataframe, remove all the bed
    # files from tmpfs to prevent a run on downloadmoreram.com
    cleanup()

    hgroup = config.feature_meta.homopolymers

    length_col = hgroup.fmt_name_len(base)
    frac_col = hgroup.fmt_name_imp_frac(base)

    merged[length_col] = merged[bed_end] - merged[bed_start] - SLOP * 2
    merged[frac_col] = 1 - (merged[PFCT_LEN_COL] / merged[length_col])
    return merged.drop(columns=[PFCT_LEN_COL])


def main(smk, config: cfg.StratoMod) -> None:
    bed_cols = config.feature_meta.bed_index
    # ASSUME this file is already sorted
    simreps = read_input(smk.input["bed"][0], bed_cols)
    merged = merge_base(
        config,
        simreps,
        cfg.Base(smk.wildcards["base"]),
        smk.input["genome"][0],
        bed_cols,
    )
    write_tsv(smk.output[0], merged, header=True)


main(snakemake, snakemake.config)  # type: ignore
