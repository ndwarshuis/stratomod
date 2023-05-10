import logging
import csv
from pathlib import Path
from typing import TypeVar
from itertools import product
from more_itertools import unzip
import pandas as pd
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from Bio import bgzf  # type: ignore


T = TypeVar("T")


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        try:
            next(bgzf.BgzfBlocks(f), None)
            return True
        except ValueError:
            return False


def read_bed(
    path: Path,
    b: cfg.BedFileParams = cfg.BedFileParams(),
    more: dict[int, cfg.PandasColumn] = {},
    chr_indices: list[cfg.ChrIndex] = [],
) -> pd.DataFrame:
    """Read a bed file as a pandas dataframe.

    Return a dataframe where the first three columns are numbered 0, 1, 2 and
    typed str, int, int (first is str regardless of how the chr names are
    formated). Columns from 'more' are appended to the end of the dataframe
    in the order given starting from 3.
    """
    bedcols = {**b.bed_cols.indexed, **more}
    df = pd.read_table(
        path,
        header=None,
        usecols=[*bedcols],
        sep=b.sep,
        comment="#",
        skiprows=b.skip_lines,
        # satisfy type checker :/
        dtype={k: v for k, v in b.bed_cols.typed.items()},
    )
    df.columns = pd.Index(bedcols.values())
    if len(chr_indices) > 0:
        f = cfg.ChrFilter(b.chr_prefix, chr_indices)
        return filter_sort_bed(f, df)
    else:
        return df


def write_bed(path: Path, df: pd.DataFrame) -> None:
    """Write a bed file in bgzip format from a dataframe.

    Dataframe is not checked to make sure it is a "real" bed file.
    """
    with bgzf.open(path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        for r in df.itertuples(index=False):
            w.writerow(r)


def filter_sort_bed(cfilt: cfg.ChrFilter, df: pd.DataFrame, n: int = 3) -> pd.DataFrame:
    from_map = {i.chr_name_full(cfilt.prefix): i.value for i in cfilt.indices}
    return filter_sort_bed_inner(from_map, df, n)


def filter_sort_bed_inner(
    from_map: dict[str, int],
    df: pd.DataFrame,
    n: int = 3,
) -> pd.DataFrame:
    chr_col = df.columns.tolist()[0]
    df[chr_col] = df[chr_col].map(from_map)
    df = df.dropna(subset=[chr_col]).astype({chr_col: int})
    return sort_bed_numerically(df, n)


def filter_chromosomes(
    chr_indices: list[cfg.ChrIndex], df: pd.DataFrame
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


def sort_bed_numerically(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that all are
    integer typed. Use 'n = 2' to sort only by chr/start, and 'n=1' to sort only
    by chr.

    """
    cols = df.columns.tolist()
    bycols = [cols[i] for i in range(0, n)]
    return df.sort_values(
        by=bycols,
        axis=0,
        ignore_index=True,
    )


# def sort_bed_numerically(df: pd.DataFrame, drop_chr: bool = True) -> pd.DataFrame:
#     # ASSUME: the first three columns correspond to a bed file and the first
#     # column has already been standardized (eg all chromosomes are numbered 1 to
#     # 24 and there are no incomplete chromosomes)
#     cols = df.columns.tolist()

#     logging.info("Numerically sorting bed")
#     return df.sort_values(
#         by=[cols[0], cols[1], cols[2]],
#         axis=0,
#         ignore_index=True,
#     )


# def read_bed_df(
#     path: str,
#     bed_mapping: dict[int, cfg.PandasColumn],
#     col_mapping: dict[int, cfg.PandasColumn],
#     chr_filter: cfg.ChrFilter,
# ) -> pd.DataFrame:
#     mapping = {**bed_mapping, **col_mapping}
#     df = pd.read_table(
#         path,
#         header=None,
#         usecols=[*mapping],
#         names=[*mapping.values()],
#     )
#     # [[*mapping]].rename(columns=mapping)
#     chr_col = df.columns.tolist()[0]
#     df_standardized = standardize_chr_column(
#         chr_filter.prefix,
#         chr_col,
#         df.astype({chr_col: str}),
#     )
#     return sort_bed_numerically(
#         filter_chromosomes(
#             chr_filter.indices,
#             df_standardized,
#         )
#     )


def merge_and_apply_stats(
    fconf: cfg.MergedFeatureGroup[T],
    bed_df: pd.DataFrame,
) -> tuple[bt, list[str]]:
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
    full_headers = [*cfg.BED_COLS, fconf.fmt_count_feature(), *headers]

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
