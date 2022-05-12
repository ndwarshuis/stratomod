import json
import pandas as pd
import numpy as np
from functools import partial
from more_itertools import duplicates_everseen
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.config import fmt_vcf_feature

logger = setup_logging(snakemake.log[0])

_fmt_vcf_feature = partial(fmt_vcf_feature, snakemake.config)

TP_LABEL = "tp"
TN_LABEL = "tn"
FP_LABEL = "fp"
FN_LABEL = "fn"

FILTERED_VAL = "RefCall"


def read_inputs(paths, input_col):
    eps = [*enumerate(paths)]
    return (
        pd.concat([read_tsv(p).assign(**{input_col: i}) for i, p in eps]),
        {p: i for i, p in eps},
    )


def process_series(opts, ser):
    trans = opts["transform"]
    _ser = pd.to_numeric(ser, errors="coerce")
    if trans == "binary":
        return (~_ser.isnull()).astype(int)
    else:
        fillval = opts["fill_na"]
        return (np.log(_ser) if trans == "log" else _ser).fillna(fillval)


def process_chr(ser):
    return (
        ser.replace({"chrX": "chr23", "chrY": "chr24"})
        .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
        .astype(int)
    )


def check_columns(wanted_cols, df_cols):
    def assert_dups(xs, msg):
        dups = [*duplicates_everseen(xs)]
        assert 0 == len(dups), f"{msg}: {dups}"
        return set(xs)

    # ASSUME we already check for duplicate feature columns when the config
    # is validated
    wanted_set = set(wanted_cols)
    df_set = assert_dups(df_cols, "input dataframe has duplicate columns")

    assert (
        df_set >= wanted_set
    ), "configuration features must be a subset of columns in input dataframe"


def select_columns(features, label_col, df):
    wanted_cols = [*features, label_col]
    check_columns(wanted_cols, df.columns.tolist())
    to_rename = {k: n for k, v in features.items() if (n := v["alt_name"]) is not None}
    return df[wanted_cols].rename(columns=to_rename)


def mask_labels(filtered_are_candidates, label_col, filter_col, df):
    # if we don't want to include filtered labels (from the perspective of
    # the truth set) they all become false negatives
    def mask(row):
        if row[filter_col] == FILTERED_VAL:
            if row[label_col] == FP_LABEL:
                return TN_LABEL
            elif row[label_col] == TP_LABEL:
                return FN_LABEL
            else:
                return row[label_col]
        else:
            return row[label_col]

    if filtered_are_candidates is False:
        # use convoluted apply to avoid slicing warnings
        df[label_col] = df.apply(mask, axis=1)
    return df


def collapse_labels(error_labels, label_col, df):
    all_labels = [*error_labels, TP_LABEL]
    return df[df[label_col].apply(lambda x: x in all_labels)].assign(
        **{label_col: lambda x: (x[label_col] == TP_LABEL).astype(int)}
    )


def process_data(
    features,
    error_labels,
    filtered_are_candidates,
    chrom_col,
    filter_col,
    label_col,
    df,
):
    for col, opts in features.items():
        if col == chrom_col:
            df[col] = process_chr(df[col])
        else:
            df[col] = process_series(opts, df[col])
    # select columns after transforms to avoid pandas asking me to make a
    # deep copy (which will happen on a slice of a slice)
    return collapse_labels(
        error_labels,
        label_col,
        select_columns(
            features,
            label_col,
            mask_labels(filtered_are_candidates, label_col, filter_col, df),
        ),
    )


def main():
    ps = snakemake.params.config
    fconf = snakemake.config["features"]
    label_col = fconf["label"]
    chrom_col = fconf["index"]["chr"]
    filter_col = _fmt_vcf_feature("filter")
    input_col = _fmt_vcf_feature("input")
    raw_df, mapped_paths = read_inputs(snakemake.input, input_col)
    processed = process_data(
        ps["features"],
        ps["error_labels"],
        ps["filtered_are_candidates"],
        chrom_col,
        filter_col,
        label_col,
        raw_df,
    )
    with open(snakemake.output["paths"], "w") as f:
        json.dump(mapped_paths, f)
    write_tsv(snakemake.output["df"], processed)


main()
