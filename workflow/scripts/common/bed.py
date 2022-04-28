import logging
from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt

BED_CHR = "chrom"
BED_START = "chromStart"
BED_END = "chromEnd"

# NOTE count and count_distinct are not allowed since those will be the same
# in multiple columns and thus are better handled using a different code path
ALLOWED_STATS = [
    "sum",
    "min",
    "max",
    "absmin",
    "absmax",
    "mean",
    "median",
    "collapse",
    "distinct",
]


def filter_chromosomes(df, chr_filter):
    if len(chr_filter) > 0:
        logging.info("Pre-filtering chromosomes: %s", ", ".join(chr_filter))
        return df[df.iloc[:, 0].isin(chr_filter)]
    return df


def sort_bed_numerically(df, drop_chr=True):
    cols = df.columns.tolist()
    tmp = "tmp_n"

    def log_unique(msg, df):
        logging.info("%s: %s", msg, ", ".join(df[cols[0]].unique().tolist()))

    def log_nrows(msg, df):
        logging.info("%s: %s", msg, df.shape[0])

    df[tmp] = (
        df[cols[0]]
        .replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
    )
    if drop_chr is True:
        logging.info("Filtering bed for complete chomosomes")
        log_nrows("Number of entries before filtering", df)
        log_unique("Unique chromosomes before filtering", df)
        df = df.dropna(axis=0, subset=[tmp])
        log_nrows("Number of entries before filtering", df)
        log_unique("Unique chromosomes after filtering", df)

    logging.info("Numerically sorting bed")
    df = (
        df.astype({tmp: int})
        .sort_values(by=[tmp, cols[1], cols[2]], axis=0, ignore_index=True)
        .drop(columns=[tmp])
    )

    return df


def merge_and_apply_stats(merge_stats, prefix, bed_df):
    assert set(merge_stats) <= set(ALLOWED_STATS), "Bad stat found"

    # compute stats on all columns except the first 3
    drop_n = 3
    stat_cols = bed_df.columns.tolist()[drop_n:]

    logging.info("Computing stats for columns: %s\n", ", ".join(stat_cols))
    logging.info("Stats to compute: %s\n", ", ".join(merge_stats))

    cols, opts, headers = unzip(
        (i + drop_n + 1, m, f"{s}_{m}")
        for (i, s), m in product(enumerate(stat_cols), merge_stats)
    )

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ["count", *opts]
    full_cols = [drop_n + 1, *cols]
    full_headers = [BED_CHR, BED_START, BED_END, f"{prefix}_count", *headers]

    logging.info("Merging regions")
    # TODO there might be a way to make pybedtools echo what it is doing, but
    # for now this is a sanity check that this crazy command is executed
    # correctly
    logging.info(
        "Using command: 'bedtools merge -i <file> -c %s -o %s'",
        ", ".join(map(str, full_cols)),
        ", ".join(full_opts),
    )

    return (
        bt.from_dataframe(bed_df).merge(c=full_cols, o=full_opts),
        full_headers,
    )
