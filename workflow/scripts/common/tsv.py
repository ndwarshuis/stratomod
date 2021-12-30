import sys
import pandas as pd


def read_tsv_from(path):
    return pd.read_csv(path, sep="\t")


def read_tsv(path):
    if path is None:
        return read_tsv_from(sys.stdin)
    # with open(path, "rt") as f:
    return read_tsv_from(path)


def write_tsv_to(f, df):
    df.to_csv(f, sep="\t", index=False, header=True)


def write_tsv(path, df):
    if path is None:
        return write_tsv_to(sys.stdout, df)
    # with open(path, "wt") as f:
    return write_tsv_to(path, df)
