from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically, standardize_chr_column
from common.cli import setup_logging
import common.config as cfg


logger = setup_logging(snakemake.log[0])


def main() -> None:
    prefix = cfg.refsetkey_to_chr_prefix(
        snakemake.config,
        snakemake.wildcards["refset_key"],
    )
    df = read_tsv(snakemake.input[0], comment="#", header=None)
    standardized = standardize_chr_column(prefix, 0, df)
    write_tsv(
        snakemake.output[0],
        sort_bed_numerically(standardized),
        header=False,
    )


main()
