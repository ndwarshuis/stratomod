import logging
import csv
import gzip
from pathlib import Path
from typing import TypeVar, Callable, TextIO
from itertools import product
from more_itertools import unzip
import pandas as pd
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from Bio import bgzf  # type: ignore


X = TypeVar("X")


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        try:
            next(bgzf.BgzfBlocks(f), None)
            return True
        except ValueError:
            return False


def with_bgzip_maybe(f: Callable[[TextIO, TextIO], X], i: str, o: str) -> X:
    # bgzf only understands latin1, so read everything as such
    hi = (
        gzip.open(i, "rt", encoding="latin1")
        if i.endswith(".gz")
        else open(i, "rt", encoding="latin1")
    )
    ho = bgzf.open(o, "wt") if o.endswith(".gz") else open(o, "wt")
    with hi as fi, ho as fo:
        return f(fi, fo)


def read_bed(
    path: Path,
    b: cfg.BedFileParams = cfg.BedFileParams(),
    more: dict[int, str] = {},
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


def sort_bed_numerically(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that all are
    integer typed. Use 'n = 2' to sort only by chr/start, and 'n=1' to sort only
    by chr.

    """
    cols = df.columns.tolist()
    bycols = [cols[i] for i in range(0, n)]
    return df.sort_values(by=bycols, axis=0, ignore_index=True)


def merge_and_apply_stats(
    fconf: cfg.MergedFeatureGroup[X],
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
    full_headers: list[str] = [*cfg.BED_COLS, fconf.count_feature[0], *headers]

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
