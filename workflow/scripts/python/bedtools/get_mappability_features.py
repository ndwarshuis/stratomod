import pandas as pd
import common.config as cfg
from typing import Any
from functools import partial
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import read_tsv, write_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])  # type: ignore


# TODO apply chr filter to these
def read_plain_bed(
    config: cfg.StratoMod,
    bedconf: cfg.BedFile,
    path: str,
    col: str,
) -> pd.DataFrame:
    logger.info("Reading mappability feature: %s", col)
    # these are just plain bed files with no extra columns
    bed_index = config.feature_names.bed_index
    bed_cols = bed_index.bed_cols_ordered()
    df = read_tsv(path, comment="#", names=bed_cols)
    # add a new column with all '1' (this will be a binary feature)
    df[col] = 1
    return standardize_chr_column(bedconf.chr_prefix, bed_index.chr, df)


def main(smk: Any, config: cfg.StratoMod) -> None:
    mapconf = config.refsetkey_to_ref(
        cfg.RefsetKey(smk.wildcards["refset_key"])
    ).annotations.mappability
    mapmeta = config.feature_names.mappability

    read_bed = partial(read_plain_bed, config)

    high = read_bed(mapconf.high, smk.input["high"][0], mapmeta.high)
    low = read_bed(mapconf.low, smk.input["low"][0], mapmeta.low)
    # subtract high from low (since the former is a subset of the latter)
    new_low = (
        bt.from_dataframe(low)
        .subtract(bt.from_dataframe(high))
        .to_dataframe(names=low.columns.tolist())
    )
    write_tsv(smk.output["high"], high)
    write_tsv(smk.output["low"], new_low)


main(snakemake, snakemake.config)  # type: ignore
