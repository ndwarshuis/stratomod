from pybedtools import BedTool as bt
from common.tsv import write_tsv
from common.bed import read_bed_df, merge_and_apply_stats
from common.cli import setup_logging
from common.config import lookup_bed_cols, bed_cols_indexed

# This database is documented here:
# http://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg38&hgta_doSchemaTable=genomicSuperDups

# ASSUME segdups dataframe is fed into this script with the chromosome column
# standardized. The column numbers below are dictionary values, and the
# corresponding feature names are the dictionary keys. Note that many feature
# names don't match the original column names in the database.

logger = setup_logging(snakemake.log[0])


def read_segdups(path, fconf, bed_cols):
    cols = fconf["columns"]
    feature_cols = {
        18: cols["alignL"],
        27: cols["fracMatchIndel"],
    }
    bed_mapping = bed_cols_indexed([1, 2, 3], bed_cols)
    return read_bed_df(path, bed_mapping, feature_cols, snakemake.params["filt"])


def merge_segdups(df, fconf, bed_cols):
    bed, names = merge_and_apply_stats(
        fconf["operations"],
        bed_cols,
        fconf["prefix"],
        df,
    )
    return bed.to_dataframe(names=names)


def main():
    fconf = snakemake.config["features"]["segdups"]
    bed_cols = lookup_bed_cols(snakemake.config)
    repeat_df = read_segdups(snakemake.input[0], fconf, bed_cols)
    merged_df = merge_segdups(repeat_df, fconf, bed_cols)
    write_tsv(snakemake.output[0], merged_df, header=True)


main()
