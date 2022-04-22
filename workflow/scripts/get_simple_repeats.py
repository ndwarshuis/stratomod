from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically, filter_chromosomes
from common.cli import setup_logging

logger = setup_logging(snakemake.log[0])

# columns from here: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
# ASSUME the input dataframe is the table downloaded as-is from the above src

BED_CHR = "chrom"
BED_START = "chromStart"
BED_END = "chromEnd"

BED_COLUMNS = {
    BED_CHR: 1,
    BED_START: 2,
    BED_END: 3,
}

FEATURE_COLUMNS = {
    "period": 5,
    "copyNum": 6,
    "consensusSize": 7,
    "perMatch": 8,
    "perIndel": 9,
    "score": 10,
    "A": 11,
    "C": 12,
    "G": 13,
    "T": 14,
}

ALL_COLUMNS = {**BED_COLUMNS, **FEATURE_COLUMNS}

MERGE_STATS = ["max", "min", "median"]
LEN_FEATURE = "region_length"


def read_simple_repeats(path):
    df = read_tsv(path, header=None)[[*ALL_COLUMNS.values()]]
    df.columns = [*ALL_COLUMNS]
    df = filter_chromosomes(df, BED_CHR, snakemake.params["filt"])
    df["AT"] = df["A"] + df["T"]
    df["GC"] = df["G"] + df["C"]
    return sort_bed_numerically(df)


def merge_simple_repeats(gfile, df):
    # compute stats on all columns except the first 3
    drop_n = 3
    stat_cols = df.columns.tolist()[drop_n:]

    logger.info("Computing stats for columns: %s\n", ", ".join(stat_cols))
    logger.info("Stats to compute: %s\n", ", ".join(MERGE_STATS))

    cols, opts, headers = unzip(
        (i + drop_n + 1, m, f"{s}_{m}")
        for (i, s), m in product(enumerate(stat_cols), MERGE_STATS)
    )

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ["count", *opts]
    full_cols = [drop_n + 1, *cols]
    full_headers = [*BED_COLUMNS, "count", *headers]

    logger.info("Merging repeat regions.")
    # TODO there might be a way to make pybedtools echo what it is doing, but
    # for now this is a sanity check that this crazy command is executed
    # correctly
    logger.info(
        "Using command: 'bedtools merge -i <file> -c %s -o %s'",
        ", ".join(map(str, full_cols)),
        ", ".join(full_opts),
    )

    merged_df = (
        bt.from_dataframe(df)
        .merge(c=full_cols, o=full_opts)
        .slop(b=5, g=gfile)
        .to_dataframe(names=full_headers)
    )
    # use the original dataframe to get the region length since we added slop to
    # the merged version
    merged_df[LEN_FEATURE] = df[BED_END] - df[BED_START]

    return merged_df


def main():
    i = snakemake.input
    repeat_df = read_simple_repeats(i.src[0])
    merged_repeat_df = merge_simple_repeats(i.genome[0], repeat_df)
    write_tsv(snakemake.output[0], merged_repeat_df, header=True)


main()
