import pandas as pd
import common.config as cfg
from common.tsv import write_tsv, read_tsv
from common.bed import standardize_chr_column
from common.cli import setup_logging


logger = setup_logging(snakemake.log[0])  # type: ignore

CHR_COL = 0


# NOTE only the first two columns matter here (chrom and length)
def read_input(path: str) -> pd.DataFrame:
    return read_tsv(path, header=None)[[0, 1]]


def sort_genome(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(by=[str(CHR_COL)], axis=0, ignore_index=True)


def write_output(path: str, df: pd.DataFrame) -> None:
    write_tsv(path, df, header=None)


def main(smk, config: cfg.StratoMod) -> None:
    prefix = cfg.refsetkey_to_ref(config, smk.wildcards["refset_key"]).genome.chr_prefix
    df = sort_genome(standardize_chr_column(prefix, CHR_COL, read_input(smk.input[0])))
    write_tsv(smk.output[0], df, header=None)


main(snakemake, snakemake.config)  # type: ignore
