import sys
import pandas as pd


def read_tsv_from(path, *args, **kwargs):
    return pd.read_csv(path, sep="\t", *args, **kwargs)


def read_tsv(path, *args, **kwargs):
    return read_tsv_from(sys.stdin if path is None else path, *args, **kwargs)


def write_tsv_to(f, df, *args, index=False, **kwargs):
    df.to_csv(f, sep="\t", index=index, *args, **kwargs)


def write_tsv(path, df, *args, **kwargs):
    return write_tsv_to(sys.stdout if path is None else path, df, *args, **kwargs)
