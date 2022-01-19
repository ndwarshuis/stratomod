import pandas as pd
from common.cli import printerr


def print_unique(msg, col, df):
    printerr("%s: %s" % (msg, ", ".join(df[col].unique().tolist())))


def print_nrows(msg, df):
    printerr("%s: %i" % (msg, df.shape[0]))


def sort_bed_numerically(df):
    cols = df.columns.tolist()
    tmp = "tmp_n"

    printerr("Filtering bed for complete chomosomes and numerically sorting")
    print_nrows("Number of entries before filtering", df)
    print_unique("Unique chromosomes before filtering", cols[0], df)

    df[tmp] = (
        df[cols[0]]
        .replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
    )
    df = df.dropna(axis=0, subset=[tmp]).astype({tmp: int})
    df = df.sort_values(by=[tmp, cols[1], cols[2]], axis=0, ignore_index=True).drop(
        columns=[tmp]
    )

    print_nrows("Number of entries after filtering", df)
    print_unique("Unique chromosomes after filtering", cols[0], df)
    print("")

    return df
