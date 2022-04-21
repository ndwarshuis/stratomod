import pandas as pd
from sys import stderr, stdout


def print_unique(msg, col, df, dev):
    unique = ", ".join(df[col].unique().tolist())
    print(f"{msg}: {unique}", file=dev)


def print_nrows(msg, df, dev):
    print(f"{msg}: {df.shape[0]}", file=dev)


def sort_bed_numerically(df, filter_chr=True, print_stderr=False):
    cols = df.columns.tolist()
    tmp = "tmp_n"
    dev = stderr if print_stderr else stdout

    df[tmp] = (
        df[cols[0]]
        .replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
    )
    if filter_chr is True:
        print("Filtering bed for complete chomosomes", file=dev)
        print_nrows("Number of entries before filtering", df, dev)
        print_unique("Unique chromosomes before filtering", cols[0], df, dev)
        df = df.dropna(axis=0, subset=[tmp])
        print_unique("Unique chromosomes after filtering", cols[0], df, dev)
        print_nrows("Number of entries before filtering", df, dev)

    print("Numerically sorting bed", file=dev)
    df = (
        df.astype({tmp: int})
        .sort_values(by=[tmp, cols[1], cols[2]], axis=0, ignore_index=True)
        .drop(columns=[tmp])
    )

    print("", file=dev)

    return df
