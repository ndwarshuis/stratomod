import logging
import pandas as pd
from common.tsv import read_tsv
from itertools import product
from more_itertools import unzip
from pybedtools import BedTool as bt
from common.config import (
    fmt_count_feature,
    fmt_merged_feature,
    bed_cols_ordered,
    fmt_strs,
)

# BED_CHR = "chrom"
# BED_START = "chromStart"
# BED_END = "chromEnd"


def filter_chromosomes(df, chr_filter):
    if len(chr_filter) > 0:
        logging.info(
            "Pre-filtering chromosomes: %s",
            fmt_strs(map(str, chr_filter)),
        )
        return df[df.iloc[:, 0].isin(chr_filter)].copy()
    return df


def standardize_chr_series(ser):
    _ser = ser.str.replace("chr", "")
    _ser[_ser == "X"] = "23"
    _ser[_ser == "Y"] = "24"
    return pd.to_numeric(_ser, errors="coerce").astype("Int64")


def sort_bed_numerically(df, drop_chr=True):
    # ASSUME: the first three columns correspond to a bed file and the first
    # column has already been standardized (eg all chromosomes are numbered 1 to
    # 24 and there are no incomplete chromosomes)
    cols = df.columns.tolist()

    # def log_unique(msg, df):
    #     logging.info("%s: %s", msg, ", ".join(df[cols[0]].unique().tolist()))

    # def log_nrows(msg, df):
    #     logging.info("%s: %s", msg, df.shape[0])

    # if drop_chr is True:
    #     logging.info("Filtering bed for complete chomosomes")
    #     log_nrows("Number of entries before filtering", df)
    #     log_unique("Unique chromosomes before filtering", df)
    #     df = df.dropna(axis=0, subset=[cols[0]])
    #     log_nrows("Number of entries before filtering", df)
    #     log_unique("Unique chromosomes after filtering", df)

    logging.info("Numerically sorting bed")
    return df.sort_values(
        by=[cols[0], cols[1], cols[2]],
        axis=0,
        ignore_index=True,
    )


def read_bed_df(path, bed_mapping, col_mapping, filt):
    mapping = {**bed_mapping, **col_mapping}
    dtypes = {0: int, 1: int, 2: int}
    df = read_tsv(path, header=None, dtype=dtypes)[[*mapping]].rename(columns=mapping)
    return sort_bed_numerically(filter_chromosomes(df, filt))


def merge_and_apply_stats(merge_stats, bed_cols, prefix, bed_df):
    # compute stats on all columns except the first 3
    drop_n = 3
    stat_cols = bed_df.columns.tolist()[drop_n:]

    logging.info("Computing stats for columns: %s\n", ", ".join(stat_cols))
    logging.info("Stats to compute: %s\n", ", ".join(merge_stats))

    cols, opts, headers = unzip(
        (i + drop_n + 1, m, fmt_merged_feature(prefix, s, m))
        for (i, s), m in product(enumerate(stat_cols), merge_stats)
    )

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ["count", *opts]
    full_cols = [drop_n + 1, *cols]
    full_headers = [
        *bed_cols_ordered(bed_cols),
        fmt_count_feature(prefix),
        *headers,
    ]

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
