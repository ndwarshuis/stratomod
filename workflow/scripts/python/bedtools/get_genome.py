import pandas as pd
import common.config as cfg
from common.tsv import write_tsv, read_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])  # type: ignore

COLS = {0: "chr", 1: "len"}


# NOTE only the first two columns matter here (chrom and length)
def read_input(path: str) -> pd.DataFrame:
    return read_tsv(path, header=None)[[*COLS]].rename(columns=COLS)


def sort_genome(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(by=[COLS[0]], axis=0, ignore_index=True)


def write_output(path: str, df: pd.DataFrame) -> None:
    write_tsv(path, df, header=None)


def main(smk, config: cfg.StratoMod) -> None:
    prefix = config.refsetkey_to_ref(smk.wildcards["refset_key"]).genome.chr_prefix
    df = sort_genome(standardize_chr_column(prefix, COLS[0], read_input(smk.input[0])))
    write_tsv(smk.output[0], df, header=None)


main(snakemake, snakemake.config)  # type: ignore
