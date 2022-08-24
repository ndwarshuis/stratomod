import logging
from functools import partial
from itertools import product
from more_itertools import unzip
import pandas as pd
from common.tsv import read_tsv
from common.functional import compose
from common.config import (
    fmt_count_feature,
    fmt_merged_feature,
    bed_cols_ordered,
    fmt_strs,
)


def filter_chromosomes(chr_filter, df):
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


def standardize_chr_column(chr_col, df):
    logging.info("Standardizing chromosome column: %s", chr_col)
    df[chr_col] = standardize_chr_series(df[chr_col])
    logging.info(
        "Removing %i rows with non-standard chromosomes",
        df[chr_col].isna().sum(),
    )
    return df.dropna(subset=[chr_col]).astype({chr_col: int})


def sort_bed_numerically(df, drop_chr=True):
    # ASSUME: the first three columns correspond to a bed file and the first
    # column has already been standardized (eg all chromosomes are numbered 1 to
    # 24 and there are no incomplete chromosomes)
    cols = df.columns.tolist()

    logging.info("Numerically sorting bed")
    return df.sort_values(
        by=[cols[0], cols[1], cols[2]],
        axis=0,
        ignore_index=True,
    )


def read_bed_df(path, bed_mapping, col_mapping, filt):
    mapping = {**bed_mapping, **col_mapping}
    df = read_tsv(path, header=None)[[*mapping]].rename(columns=mapping)
    chr_col = df.columns.tolist()[0]
    return compose(
        sort_bed_numerically,
        partial(filter_chromosomes, filt),
        partial(standardize_chr_column, chr_col),
    )(df.astype({chr_col: str}))


def merge_and_apply_stats(merge_stats, bed_cols, prefix, bed_df):
    # import this here so we can import other functions in this module
    # without pulling in bedtools
    from pybedtools import BedTool as bt

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
