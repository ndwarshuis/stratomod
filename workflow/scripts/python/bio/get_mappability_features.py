from pathlib import Path
import pandas as pd
import common.config as cfg
from typing import Any
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import write_tsv
from common.bed import read_bed
from common.io import setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, config: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    cs = config.refsetkey_to_chr_indices(rsk)
    mapconf = config.refsetkey_to_ref(rsk).annotations.mappability
    mapmeta = config.feature_names.mappability

    def read_map_bed(p: Path, ps: cfg.BedFileParams, col: str) -> pd.DataFrame:
        logger.info("Reading mappability feature: %s", col)
        df = read_bed(p, ps, {}, cs)
        df[col] = 1
        return df

    high = read_map_bed(smk.input["high"][0], mapconf.high.params, mapmeta.high)
    low = read_map_bed(smk.input["low"][0], mapconf.low.params, mapmeta.low)
    # subtract high from low (since the former is a subset of the latter)
    new_low = (
        bt.from_dataframe(low)
        .subtract(bt.from_dataframe(high))
        .to_dataframe(names=low.columns.tolist())
    )
    write_tsv(smk.output["high"], high)
    write_tsv(smk.output["low"], new_low)


main(snakemake, snakemake.config)  # type: ignore
