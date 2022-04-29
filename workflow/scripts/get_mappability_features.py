from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv

# the only thing this script does (so far) is subtract the high mappability
# regions from the low regions (since high is a subset of low)


def main():
    high = bt.from_dataframe(read_tsv(snakemake.input["high"][0]))
    low = read_tsv(snakemake.input["low"][0])
    names = low.columns.tolist()
    new = bt.from_dataframe(low).subtract(high).to_dataframe(names=names)
    write_tsv(snakemake.output[0], new)


main()
