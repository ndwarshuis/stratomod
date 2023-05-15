from pathlib import Path
import pandas as pd
import common.config as cfg
from typing import Any, cast
from pybedtools import BedTool as bt  # type: ignore
from pybedtools import cleanup
from common.tsv import write_tsv
from common.io import setup_logging
from common.bed import read_bed


logger = setup_logging(snakemake.log[0])  # type: ignore

# temporary columns used for dataframe processing
BASE_COL = "_base"
PFCT_LEN_COL = "_perfect_length"

SLOP = 1


def read_input(path: Path) -> pd.DataFrame:
    logger.info("Reading dataframe from %s", path)
    return read_bed(path, more={3: cfg.PandasColumn(BASE_COL)})


def merge_base(
    config: cfg.StratoMod,
    df: pd.DataFrame,
    base: cfg.Base,
    genome: str,
) -> pd.DataFrame:
    logger.info("Filtering bed file for %ss", base)
    _df = df[df[BASE_COL] == f"unit={base.value}"].drop(columns=[BASE_COL])
    logger.info("Merging %s rows for %ss", len(_df), base)
    # Calculate the length of each "pure" homopolymer (eg just "AAAAAAAA").
    # Note that this is summed in the merge below, and the final length based
    # on start/end won't necessarily be this sum because of the -d 1 parameter
    _df[PFCT_LEN_COL] = _df[cfg.BED_END] - _df[cfg.BED_START]
    merged = cast(
        pd.DataFrame,
        bt.from_dataframe(_df)
        .merge(d=1, c=[4], o=["sum"])
        .slop(b=SLOP, g=genome)
        .to_dataframe(names=[*cfg.BED_COLS, PFCT_LEN_COL]),
    )
    # these files are huge; now that we have a dataframe, remove all the bed
    # files from tmpfs to prevent a run on downloadmoreram.com
    cleanup()

    hgroup = config.feature_definitions.homopolymers

    length_col = hgroup.fmt_name(base, lambda x: x.len)
    frac_col = hgroup.fmt_name(base, lambda x: x.imp_frac)

    merged[length_col] = merged[cfg.BED_END] - merged[cfg.BED_START] - SLOP * 2
    merged[frac_col] = 1 - (merged[PFCT_LEN_COL] / merged[length_col])
    return merged.drop(columns=[PFCT_LEN_COL])


def main(smk: Any, config: cfg.StratoMod) -> None:
    # ASSUME this file is already sorted
    simreps = read_input(smk.input["bed"][0])
    merged = merge_base(
        config,
        simreps,
        cfg.Base(smk.wildcards["base"]),
        smk.input["genome"][0],
    )
    write_tsv(smk.output[0], merged, header=True)


main(snakemake, snakemake.config)  # type: ignore
