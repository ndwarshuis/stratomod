import sys
import pandas as pd
from typing import Any, cast


def read_tsv_from(path: str, *args: Any, **kwargs: Any) -> pd.DataFrame:
    return cast(pd.DataFrame, pd.read_csv(path, sep="\t", *args, **kwargs))


def read_tsv(path: str, *args: Any, **kwargs: Any) -> pd.DataFrame:
    return read_tsv_from(sys.stdin if path is None else path, *args, **kwargs)


def write_tsv_to(
    f: str,
    df: pd.DataFrame,
    *args: Any,
    index: bool = False,
    **kwargs: Any,
) -> None:
    df.to_csv(f, sep="\t", index=index, *args, **kwargs)  # type: ignore


def write_tsv(path: str, df: pd.DataFrame, *args: Any, **kwargs: Any) -> None:
    return write_tsv_to(sys.stdout if path is None else path, df, *args, **kwargs)
