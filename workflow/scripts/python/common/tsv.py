import pandas as pd


def write_tsv(path: str, df: pd.DataFrame, header: bool | list[str] = True) -> None:
    df.to_csv(path, sep="\t", header=header)
