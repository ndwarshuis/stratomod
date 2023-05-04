import json
import pandas as pd
import numpy as np
from numpy.typing import NDArray
from typing import Any, Hashable, cast, TypedDict
from common.io import setup_logging
from common.ebm import read_model
from common.tsv import write_tsv
import common.config as cfg
from interpret.glassbox import ExplainableBoostingClassifier  # type: ignore
from enum import Enum

setup_logging(snakemake.log[0])  # type: ignore


IndexedVectors = dict[int, NDArray[np.float64]]
NamedVectors = dict[str, NDArray[np.float64]]

EBMUniData = TypedDict(
    "EBMUniData",
    {"type": str, "names": list[float | str], "scores": NDArray[np.float64]},
)


class VarType(Enum):
    INT = "interaction"
    CNT = "continuous"
    CAT = "categorical"


AllFeatures = dict[str, tuple[VarType, int]]


Variable = TypedDict("Variable", {"name": str, "type": str})


BivariateData = TypedDict(
    "BivariateData",
    {
        "left": Variable,
        "right": Variable,
        "df": dict[Hashable, float],
    },
)


GlobalScoreData = TypedDict(
    "GlobalScoreData",
    {
        "variable": list[str],
        "score": list[float],
    },
)


UnivariateDF = TypedDict(
    "UnivariateDF",
    {
        "value": list[str | float],
        "score": list[float],
        "stdev": list[float],
    },
)


UnivariateData = TypedDict(
    "UnivariateData",
    {
        "name": str,
        "vartype": str,
        "df": UnivariateDF,
    },
)


ModelData = TypedDict(
    "ModelData",
    {
        "global_scores": GlobalScoreData,
        "intercept": float,
        "univariate": list[UnivariateData],
        "bivariate": list[BivariateData],
    },
)


# TODO there is no reason this can't be done immediately after training
# just to avoid the pickle thing


def array_to_list(arr: NDArray[np.float64], repeat_last: bool) -> list[float]:
    # cast needed since this can return a nested list depending on number of dims
    al = cast(list[float], arr.tolist())
    return al + [al[-1]] if repeat_last else al


def get_univariate_df(
    continuous: bool,
    feature_data: EBMUniData,
    stdev: NDArray[np.float64],
) -> UnivariateDF:
    def proc_scores(scores: NDArray[np.float64]) -> list[float]:
        return array_to_list(scores, continuous)

    return UnivariateDF(
        value=feature_data["names"],
        score=proc_scores(feature_data["scores"]),
        # For some reason, the standard deviations array has an extra 0 in
        # in the front and thus is one longer than the scores array.
        stdev=proc_scores(stdev[1:]),
    )


def build_scores_array(
    arr: NDArray[np.float64],
    left_type: VarType,
    right_type: VarType,
) -> NDArray[np.float64]:
    # any continuous dimension is going to be one less than the names length,
    # so copy the last row/column to the end in these cases
    if left_type == VarType.CNT:
        arr = np.vstack((arr, arr[-1, :]))
    if right_type == VarType.CNT:
        arr = np.column_stack((arr, arr[:, -1]))
    return arr


def get_bivariate_df(
    all_features: AllFeatures,
    ebm_global: ExplainableBoostingClassifier,
    name: str,
    data_index: int,
    stdevs: IndexedVectors,
) -> BivariateData:
    def lookup_feature_type(name: str) -> VarType:
        return all_features[name][0]

    feature_data = ebm_global.data(data_index)
    # left is first dimension, right is second
    left_name, right_name = tuple(name.split(" x "))

    left_type = lookup_feature_type(left_name)
    right_type = lookup_feature_type(right_name)

    left_index = pd.Index(feature_data["left_names"], name="left_value")
    right_index = pd.Index(feature_data["right_names"], name="right_value")

    def stack_array(arr: NDArray[np.float64], name: str) -> "pd.Series[float]":
        return cast(
            "pd.Series[float]",
            pd.DataFrame(
                build_scores_array(arr, left_type, right_type),
                index=left_index,
                columns=right_index,
            ).stack(),
        ).rename(name)

    # the standard deviations are in an array that has 1 larger shape than the
    # scores array in both directions where the first row/column is all zeros.
    # Not sure why it is all zeros, but in order to make it line up with the
    # scores array we need to shave off the first row/column.
    return BivariateData(
        left=Variable(name=left_name, type=left_type.value),
        right=Variable(name=right_name, type=right_type.value),
        df=pd.concat(
            [
                stack_array(feature_data["scores"], "score"),
                stack_array(stdevs[data_index][1:, 1:], "stdev"),
            ],
            axis=1,
        )
        .reset_index()
        .to_dict(orient="list"),
    )


def get_global_scores(ebm_global: ExplainableBoostingClassifier) -> GlobalScoreData:
    glob = ebm_global.data()
    return GlobalScoreData(variable=glob["names"], score=glob["scores"])


def get_univariate_list(
    ebm_global: ExplainableBoostingClassifier,
    all_features: AllFeatures,
    stdevs: IndexedVectors,
) -> list[UnivariateData]:
    return [
        UnivariateData(
            name=name,
            vartype=vartype.value,
            df=get_univariate_df(
                vartype == VarType.CNT,
                ebm_global.data(i),
                stdevs[i],
            ),
        )
        for name, (vartype, i) in all_features.items()
        if vartype in [VarType.CNT, VarType.CAT]
    ]


def get_bivariate_list(
    ebm_global: ExplainableBoostingClassifier,
    all_features: AllFeatures,
    stdevs: IndexedVectors,
) -> list[BivariateData]:
    return [
        get_bivariate_df(all_features, ebm_global, name, i, stdevs)
        for name, (vartype, i) in all_features.items()
        if vartype == VarType.INT
    ]


def get_model(ebm: ExplainableBoostingClassifier) -> ModelData:
    ebm_global = ebm.explain_global()
    stdevs = cast(IndexedVectors, ebm.term_standard_deviations_)
    all_features = {
        cast(str, n): (VarType(t), i)
        for i, (n, t) in enumerate(
            map(tuple, ebm_global.selector[["Name", "Type"]].to_numpy())
        )
    }
    return ModelData(
        global_scores=get_global_scores(ebm_global),
        intercept=ebm.intercept_[0],
        univariate=get_univariate_list(ebm_global, all_features, stdevs),
        bivariate=get_bivariate_list(ebm_global, all_features, stdevs),
    )


def write_model_json(path: str, ebm: ExplainableBoostingClassifier) -> None:
    with open(path, "w") as f:
        json.dump(get_model(ebm), f)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    sin = smk.input
    sout = smk.output

    ebm = read_model(sin["model"])
    write_model_json(sout["model"], ebm)

    label = sconf.feature_names.label
    bed_cols = sconf.feature_names.all_index_cols()

    def write_predictions(xpath: str, ypath: str, out_path: str) -> None:
        X = pd.read_table(xpath).drop(columns=bed_cols)
        y = pd.read_table(ypath)
        y_pred = pd.DataFrame(
            {
                "prob": ebm.predict_proba(X)[::, 1],
                "label": y[label],
            }
        )
        write_tsv(out_path, y_pred)

    write_predictions(sin["train_x"], sin["train_y"], sout["train_predictions"])
    write_predictions(sin["test_x"], sin["test_y"], sout["predictions"])


main(snakemake, snakemake.config)  # type: ignore
