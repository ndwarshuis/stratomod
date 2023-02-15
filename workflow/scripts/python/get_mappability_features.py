import pandas as pd
import common.config as cfg
from functools import partial
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import read_tsv, write_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])  # type: ignore


def read_plain_bed(
    config: cfg.StratoMod,
    bedconf: cfg.BedFile,
    path: str,
    key: str,
) -> pd.DataFrame:
    logger.info("Reading mappability %s", key)
    # these are just plain bed files with no extra columns
    bed_cols = [*cfg.lookup_bed_cols(config).dict().values()]
    df = read_tsv(path, comment="#", names=bed_cols)
    # add a new column with all '1' (this will be a binary feature)
    bin_col = cfg.fmt_mappability_feature(config, key)
    df[bin_col] = 1
    return standardize_chr_column(bedconf.chr_prefix, bed_cols[0], df)


def main(smk, config: cfg.StratoMod) -> None:
    mapconf = cfg.refsetkey_to_ref(
        config, smk.wildcards["refset_key"]
    ).annotations.mappability

    read_bed = partial(read_plain_bed, config)

    high = read_bed(mapconf.high, smk.input["high"][0], "high")
    low = read_bed(mapconf.low, smk.input["low"][0], "low")
    # subtract high from low (since the former is a subset of the latter)
    new_low = (
        bt.from_dataframe(low)
        .subtract(bt.from_dataframe(high))
        .to_dataframe(names=low.columns.tolist())
    )
    write_tsv(smk.output["high"], high)
    write_tsv(smk.output["low"], new_low)


main(snakemake, snakemake.config)  # type: ignore
