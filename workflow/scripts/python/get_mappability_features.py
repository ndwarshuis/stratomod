import common.config as cfg
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])


def read_plain_bed(path, key, prefix):
    logger.info("Reading mappability %s", key)
    # these are just plain bed files with no extra columns
    bed_cols = [*cfg.lookup_bed_cols(snakemake.config).values()]
    df = read_tsv(path, comment="#", names=bed_cols)
    # add a new column with all '1' (this will be a binary feature)
    bin_col = cfg.fmt_mappability_feature(snakemake.config, key)
    df[bin_col] = 1
    return standardize_chr_column(prefix, bed_cols[0], df)


def main():
    prefix = cfg.refsetkey_to_chr_prefix(
        snakemake.config, snakemake.wildcards["refset_key"]
    )
    high = read_plain_bed(snakemake.input["high"][0], "high", prefix)
    low = read_plain_bed(snakemake.input["low"][0], "low", prefix)

    # subtract high from low (since the former is a subset of the latter)
    new_low = (
        bt.from_dataframe(low)
        .subtract(bt.from_dataframe(high))
        .to_dataframe(names=low.columns.tolist())
    )
    write_tsv(snakemake.output["high"], high)
    write_tsv(snakemake.output["low"], new_low)


main()
