import pandas as pd
import numpy as np
import common.config as cfg
from typing import Optional, cast
from typing_extensions import assert_never
from more_itertools import duplicates_everseen

# TODO don't hardcode this (and also turn into a list)
FILTERED_VAL = "RefCall"


def process_series(opts: cfg.Feature, ser: "pd.Series[float]") -> "pd.Series[float]":
    trans = opts.transform
    _ser = pd.to_numeric(ser, errors="coerce")
    if trans == cfg.Transform.BINARY:
        return (~_ser.isnull()).astype(float)
    else:
        fillval = opts.fill_na
        return cast(
            "pd.Series[float]",
            (np.log10(_ser) if trans == cfg.Transform.LOG else _ser).fillna(fillval),
        )


# def process_chr(ser):
#     return (
#         ser.replace({"chrX": "chr23", "chrY": "chr24"})
#         .str.extract(r"^chr(\d|1\d|2[0-4])$", expand=False)
#         .astype(int)
#     )


def unary_to_series(x: cfg.UnaryExpression, df: pd.DataFrame) -> "pd.Series[float]":
    f = x.function
    s = df[x.arg]
    if f is cfg.UnaryFunction.LOG:
        return cast("pd.Series[float]", np.log10(s))
    if f is cfg.UnaryFunction.BINARY:
        return (~s.isnull()).astype(float)
    else:
        assert_never(f)


def eval_equation_predicate(
    x: cfg.EquationPredicate, df: pd.DataFrame
) -> "pd.Series[bool]":
    def column(c: cfg.FeatureKey | float) -> "pd.Series[float]" | float:
        if isinstance(c, str):
            return df[c]
        elif isinstance(c, float):
            return c
        else:
            assert_never(c)

    left = column(x.left)
    right = column(x.right)

    r = x.relation
    if r is cfg.RelationalOperator.GE:
        a = np.greater_equal(left, right)
    elif r is cfg.RelationalOperator.GT:
        a = np.greater(left, right)
    elif r is cfg.RelationalOperator.EQ:
        a = np.equal(left, right)
    elif r is cfg.RelationalOperator.LE:
        a = np.less_equal(left, right)
    elif r is cfg.RelationalOperator.LT:
        a = np.less(left, right)
    elif r is cfg.RelationalOperator.NE:
        a = ~np.equal(left, right)
    else:
        assert_never(r)
    return pd.Series(a)


def eval_predicate_expression(
    x: cfg.PredicateExpression, df: pd.DataFrame
) -> "pd.Series[bool]":
    if isinstance(x, cfg.IsMissingPredicate):
        return df[x.is_missing].isnull()
    elif isinstance(x, cfg.EquationPredicate):
        return eval_equation_predicate(x, df)
    elif isinstance(x, cfg.AndPredicate):
        a = x._and
        return eval_predicate_expression(a[0], df) & eval_predicate_expression(a[1], df)
    elif isinstance(x, cfg.OrPredicate):
        o = x._or
        return eval_predicate_expression(o[0], df) | eval_predicate_expression(o[1], df)
    elif isinstance(x, cfg.NotPredicate):
        return ~eval_predicate_expression(x._not, df)
    else:
        assert_never(x)


def if_then_to_series(x: cfg.IfThenExpression, df: pd.DataFrame) -> "pd.Series[float]":
    return pd.Series(
        np.where(
            eval_predicate_expression(x.predicate, df),
            expression_to_series(x.then, df),
            expression_to_series(x._else, df),
        )
    )


def const_to_series(
    x: cfg.ConstExpression, df: pd.DataFrame
) -> "pd.Series[float]" | float:
    c = x.const
    if isinstance(c, str):
        return df[c]
    elif isinstance(c, float):
        return c
    else:
        assert_never(c)


def expression_to_series(
    x: cfg.VirtualExpression, df: pd.DataFrame
) -> "pd.Series[float]" | float:
    if isinstance(x, cfg.UnaryExpression):
        return unary_to_series(x, df)
    elif isinstance(x, cfg.IfThenExpression):
        return if_then_to_series(x, df)
    elif isinstance(x, cfg.ConstExpression):
        return const_to_series(x, df)
    else:
        assert_never(x)


def process_columns(
    features: dict[cfg.FeatureKey, cfg.Feature],
    virt_features: dict[cfg.FeatureKey, cfg.VirtualFeature],
    df: pd.DataFrame,
) -> pd.DataFrame:
    for vk, vv in virt_features.items():
        df[vk] = expression_to_series(vv.expression, df)
    for col, opts in features.items():
        df[col] = process_series(opts, df[col])
    return df


def check_columns(wanted_cols: list[cfg.FeatureKey], df_cols: list[str]) -> None:
    def assert_dups(xs: list[str], msg: str) -> set[str]:
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
    features: dict[cfg.FeatureKey, cfg.Feature],
    idx_cols: list[cfg.FeatureKey],
    label_col: Optional[cfg.FeatureKey],
    df: pd.DataFrame,
) -> pd.DataFrame:
    wanted_cols = [*features] if label_col is None else [*features, label_col]
    check_columns(wanted_cols, df.columns.tolist())
    # TODO ensure that "VCF_input" is included even when we don't "want" it,
    # since this will work with whatever the column represented by "idx_col"
    # to make a complete index mapping back to the input variant
    all_cols = idx_cols + wanted_cols
    return df[all_cols]


def mask_labels(
    filtered_are_candidates: bool,
    label_col: str,
    filter_col: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    # if we don't want to include filtered labels (from the perspective of
    # the truth set) they all become false negatives
    def mask(row: dict[str, str]) -> str:
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
    error_labels: set[cfg.ErrorLabel],
    label_col: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    all_labels = [*[x.value for x in error_labels], cfg.AnyLabel.TP.value]
    return df[df[label_col].apply(lambda x: x in all_labels)].assign(
        **{label_col: lambda x: (x[label_col] == cfg.AnyLabel.TP.value).astype(int)}
    )


def process_labeled_data(
    features: dict[cfg.FeatureKey, cfg.Feature],
    virtual_features: dict[cfg.FeatureKey, cfg.VirtualFeature],
    error_labels: set[cfg.ErrorLabel],
    filtered_are_candidates: bool,
    idx_cols: list[cfg.FeatureKey],
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
                process_columns(features, virtual_features, df),
            ),
        ),
    )


def process_unlabeled_data(
    features: dict[cfg.FeatureKey, cfg.Feature],
    virtual_features: dict[cfg.FeatureKey, cfg.VirtualFeature],
    idx_cols: list[cfg.FeatureKey],
    df: pd.DataFrame,
) -> pd.DataFrame:
    return select_columns(
        features,
        idx_cols,
        None,
        process_columns(features, virtual_features, df),
    )
