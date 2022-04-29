import yaml
import pandas as pd
import numpy as np
from more_itertools import duplicates_everseen
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging

logger = setup_logging(snakemake.log[0])


# TODO don't hardcode this
LABEL_COL = "label"
CHROM_COL = "VCF_CHROM"
FILTER_COL = "VCF_FILTER"
INPUT_COL = "VCF_input"

TP_LABEL = "tp"
TN_LABEL = "tn"
FP_LABEL = "fp"
FN_LABEL = "fn"

FILTERED_VAL = "RefCall"


def read_inputs(paths):
    eps = [*enumerate(paths)]
    return (
        pd.concat([read_tsv(p).assign(**{INPUT_COL: i}) for i, p in eps]),
        {p: i for i, p in eps},
    )


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
    wanted_cols = [*features, LABEL_COL]
    check_columns(wanted_cols, df.columns.tolist())
    return df[wanted_cols]


def mask_labels(include_filtered, df):
    # if we don't want to include filtered labels (from the perspective of
    # the truth set) they all become false negatives
    def mask(row):
        if row[FILTER_COL] == FILTERED_VAL:
            if row[LABEL_COL] == FP_LABEL:
                return TN_LABEL
            elif row[LABEL_COL] == TP_LABEL:
                return FN_LABEL
            else:
                return row[LABEL_COL]
        else:
            return row[LABEL_COL]

    if include_filtered is False:
        # use convoluted apply to avoid slicing warnings
        df[LABEL_COL] = df.apply(mask, axis=1)
    return df


def collapse_labels(error_labels, df):
    all_labels = [*error_labels, TP_LABEL]
    return df[df[LABEL_COL].apply(lambda x: x in all_labels)].assign(
        **{LABEL_COL: lambda x: (x[LABEL_COL] == TP_LABEL).astype(int)}
    )


def process_data(features, error_labels, include_filtered, df):
    for col, opts in features.items():
        if col == CHROM_COL:
            df[col] = process_chr(df[col])
        else:
            df[col] = process_series(opts, df[col])
    # select columns after transforms to avoid pandas asking me to make a
    # deep copy (which will happen on a slice of a slice)
    return collapse_labels(
        error_labels,
        select_columns(
            features,
            mask_labels(include_filtered, df),
        ),
    )


def main():
    raw_df, mapped_paths = read_inputs(snakemake.input)
    ps = snakemake.params.config
    processed = process_data(
        ps["features"],
        ps["error_labels"],
        ps["include_filtered"],
        raw_df,
    )
    with open(snakemake.output["paths"], "w") as f:
        yaml.dump(mapped_paths, f)
    write_tsv(snakemake.output["df"], processed)


main()
