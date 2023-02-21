import logging
from typing import List, Dict, Tuple, Set
from itertools import product
from more_itertools import unzip
import pandas as pd
from common.tsv import read_tsv
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore


def filter_chromosomes(
    chr_indices: Set[cfg.ChrIndex], df: pd.DataFrame
) -> pd.DataFrame:
    cs = [x.value for x in chr_indices]
    if len(cs) > 0:
        logging.info("Pre-filtering chromosomes: %s", cs)
        _df = df[df.iloc[:, 0].isin(cs)].copy()
        assert not _df.empty
        return _df
    return df


def standardize_chr_series(prefix: str, ser: "pd.Series[str]") -> "pd.Series[int]":
    _ser = ser if prefix == "" else ser.str.replace(prefix, "")
    for c in [cfg.ChrIndex.CHRX, cfg.ChrIndex.CHRY]:
        _ser[_ser == c.chr_name] = str(c.value)
    return pd.to_numeric(_ser, errors="coerce").astype("Int64")


def standardize_chr_column(
    prefix: str,
    chr_col: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    logging.info("Standardizing chromosome column: %s", chr_col)
    df[chr_col] = standardize_chr_series(prefix, df[chr_col])
    logging.info(
        "Removing %i rows with non-standard chromosomes",
        df[chr_col].isna().sum(),
    )
    return df.dropna(subset=[chr_col]).astype({chr_col: "int"})


def sort_bed_numerically(df: pd.DataFrame, drop_chr: bool = True) -> pd.DataFrame:
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


def read_bed_df(
    path: str,
    bed_mapping: Dict[int, str],
    col_mapping: Dict[int, str],
    chr_filter: cfg.ChrFilter,
) -> pd.DataFrame:
    mapping = {**bed_mapping, **col_mapping}
    df = read_tsv(path, header=None)[[*mapping]].rename(columns=mapping)
    chr_col = df.columns.tolist()[0]
    df_standardized = standardize_chr_column(
        chr_filter.prefix,
        chr_col,
        df.astype({chr_col: str}),
    )
    return sort_bed_numerically(
        filter_chromosomes(
            chr_filter.indices,
            df_standardized,
        )
    )


def merge_and_apply_stats(
    bed_cols: cfg.BedIndex,
    fconf: cfg.MergedFeatureGroup,
    bed_df: pd.DataFrame,
) -> Tuple[bt, List[str]]:
    # import this here so we can import other functions in this module
    # without pulling in bedtools

    # compute stats on all columns except the first 3
    drop_n = 3
    stat_cols = bed_df.columns.tolist()[drop_n:]

    logging.info("Computing stats for columns: %s\n", stat_cols)
    logging.info("Stats to compute: %s\n", [x.value for x in fconf.operations])

    cols, opts, headers = unzip(
        (i + drop_n + 1, m.value, fconf.fmt_merged_feature(s, m))
        for (i, s), m in product(enumerate(stat_cols), fconf.operations)
    )

    # just use one column for count since all columns will produce the same
    # number
    full_opts = ["count", *opts]
    full_cols = [drop_n + 1, *cols]
    full_headers = [
        *bed_cols.bed_cols_ordered(),
        fconf.fmt_count_feature(),
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
