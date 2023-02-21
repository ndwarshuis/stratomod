import pandas as pd
import numpy as np
import common.config as cfg
from typing import List, Set, Dict, Optional
from more_itertools import duplicates_everseen

# TODO don't hardcode this (and also turn into a list)
FILTERED_VAL = "RefCall"


def process_series(opts: cfg.Feature, ser: pd.Series[float]) -> pd.Series[float]:
    trans = opts.transform
    _ser = pd.to_numeric(ser, errors="coerce")
    if trans == "binary":
        return (~_ser.isnull()).astype(int)
    else:
        fillval = opts.fill_na
        return (np.log10(_ser) if trans == "log" else _ser).fillna(fillval)


# def process_chr(ser):
#     return (
#         ser.replace({"chrX": "chr23", "chrY": "chr24"})
#         .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
#         .astype(int)
#     )


def process_columns(
    features: Dict[cfg.FeatureKey, cfg.Feature],
    df: pd.DataFrame,
) -> pd.DataFrame:
    for col, opts in features.items():
        df[col] = process_series(opts, df[col])
    return df


def check_columns(wanted_cols: List[cfg.FeatureKey], df_cols: List[str]) -> None:
    def assert_dups(xs: List[str], msg: str) -> Set[str]:
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


def select_columns(
    features: Dict[cfg.FeatureKey, cfg.Feature],
    idx_cols: List[cfg.FeatureKey],
    label_col: Optional[cfg.FeatureKey],
    df: pd.DataFrame,
) -> pd.DataFrame:
    wanted_cols = [*features] if label_col is None else [*features, label_col]
    check_columns(wanted_cols, df.columns.tolist())
    # TODO ensure that "VCF_input" is included even when we don't "want" it,
    # since this will work with whatever the column represented by "idx_col"
    # to make a complete index mapping back to the input variant
    all_cols = idx_cols + wanted_cols
    to_rename = {k: n for k, v in features.items() if (n := v.alt_name) is not None}
    return df[all_cols].rename(columns=to_rename)


def mask_labels(
    filtered_are_candidates: bool,
    label_col: str,
    filter_col: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    # if we don't want to include filtered labels (from the perspective of
    # the truth set) they all become false negatives
    def mask(row: Dict[str, str]) -> str:
        if row[filter_col] == FILTERED_VAL:
            if row[label_col] == cfg.AnyLabel.FP.value:
                return cfg.AnyLabel.TN.value
            elif row[label_col] == cfg.AnyLabel.TP.value:
                return cfg.AnyLabel.FN.value
            else:
                return row[label_col]
        else:
            return row[label_col]

    if filtered_are_candidates is False:
        # use convoluted apply to avoid slicing warnings
        df[label_col] = df.apply(mask, axis=1)
    return df


def collapse_labels(
    error_labels: Set[cfg.ErrorLabel],
    label_col: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    all_labels = [*[x.value for x in error_labels], cfg.AnyLabel.TP.value]
    return df[df[label_col].apply(lambda x: x in all_labels)].assign(
        **{label_col: lambda x: (x[label_col] == cfg.AnyLabel.TP.value).astype(int)}
    )


def process_labeled_data(
    features: Dict[cfg.FeatureKey, cfg.Feature],
    error_labels: Set[cfg.ErrorLabel],
    filtered_are_candidates: bool,
    idx_cols: List[cfg.FeatureKey],
    filter_col: str,
    label_col: cfg.FeatureKey,
    df: pd.DataFrame,
) -> pd.DataFrame:
    # select columns after transforms to avoid pandas asking me to make a
    # deep copy (which will happen on a slice of a slice)
    return collapse_labels(
        error_labels,
        label_col,
        select_columns(
            features,
            idx_cols,
            label_col,
            mask_labels(
                filtered_are_candidates,
                label_col,
                filter_col,
                process_columns(features, df),
            ),
        ),
    )


def process_unlabeled_data(
    features: Dict[cfg.FeatureKey, cfg.Feature],
    idx_cols: List[cfg.FeatureKey],
    df: pd.DataFrame,
) -> pd.DataFrame:
    return select_columns(features, idx_cols, None, process_columns(features, df))
