import logging
import pandas as pd


def filter_chromosomes(df, col, chr_filter):
    if len(chr_filter) > 0:
        logging.info("Pre-filtering chromosomes: %s", ", ".join(chr_filter))
        return df[df[col].isin(chr_filter)]
    return df


def sort_bed_numerically(df, drop_chr=True):
    cols = df.columns.tolist()
    tmp = "tmp_n"

    def log_unique(msg, df):
        logging.info("%s: %s", msg, ", ".join(df[cols[0]].unique().tolist()))

    def log_nrows(msg, df):
        logging.info("%s: %s", msg, df.shape[0])

    df[tmp] = (
        df[cols[0]]
        .replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
    )
    if drop_chr is True:
        logging.info("Filtering bed for complete chomosomes")
        log_nrows("Number of entries before filtering", df)
        log_unique("Unique chromosomes before filtering", df)
        df = df.dropna(axis=0, subset=[tmp])
        log_nrows("Number of entries before filtering", df)
        log_unique("Unique chromosomes after filtering", df)

    logging.info("Numerically sorting bed")
    df = (
        df.astype({tmp: int})
        .sort_values(by=[tmp, cols[1], cols[2]], axis=0, ignore_index=True)
        .drop(columns=[tmp])
    )

    return df
