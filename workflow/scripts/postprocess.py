import pandas as pd
import numpy as np
from more_itertools import duplicates_everseen
from common.tsv import read_tsv, write_tsv


# TODO don't hardcode this
LABEL = "label"
TP_LABEL = "tp"
CHROM_COL = "CHROM"


def process_series(opts, ser):
    log_trans = opts["log_transform"]
    fillval = opts["fill_na"]
    _ser = pd.to_numeric(ser, errors="coerce")
    return (np.log(_ser) if log_trans else _ser).fillna(fillval)


def process_chr(ser):
    return (
        ser.replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
        .astype(int)
    )


def check_columns(wanted_cols, df_cols):
    wanted_set = set(wanted_cols)
    df_set = set(df_cols)
    assert len(wanted_set) == len(
        wanted_cols
    ), "duplicate configuration features detected: {}".format(
        list(duplicates_everseen(wanted_cols)),
    )
    assert len(df_set) == len(
        df_cols
    ), "input dataframe has duplicate columns: {}".format(
        list(duplicates_everseen(df_cols)),
    )
    assert (
        df_set >= wanted_set
    ), "configuration features must be a subset of columns in input dataframe"


def select_columns(features, df):
    wanted_cols = [*features] + [LABEL]
    check_columns(wanted_cols, df.columns.tolist())
    return df[wanted_cols]


def collapse_labels(error_labels, df):
    all_labels = [*error_labels, TP_LABEL]
    return df[df[LABEL].apply(lambda x: x in all_labels)].assign(
        **{LABEL: lambda x: (x[LABEL] == TP_LABEL).astype(int)}
    )


def process_data(features, error_labels, df):
    for col, opts in features.items():
        if col == CHROM_COL:
            df[col] = process_chr(df[col])
        else:
            df[col] = process_series(opts, df[col])
    # select columns after transforms to avoid pandas asking me to make a
    # deep copy (which will happen on a slice of a slice)
    return collapse_labels(error_labels, select_columns(features, df))


def main():
    raw = read_tsv(snakemake.input[0])
    ps = snakemake.params
    processed = process_data(ps.features, ps.error_labels, raw)
    write_tsv(snakemake.output[0], processed)


main()
