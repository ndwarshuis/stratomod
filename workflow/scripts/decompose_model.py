import pickle
import json
import pandas as pd
import numpy as np
from dash import html
from sklearn import metrics
from sklearn.model_selection import train_test_split
from interpret import set_visualize_provider
from interpret.provider import InlineProvider
from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

# TODO some of these imports probably aren't actually needed
# TODO also, there is no reason this can't be done immediately after training
# just to avoid the pickle thing


def read_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def array_to_list(arr, repeat_last):
    al = arr.tolist()
    return al + [al[-1]] if repeat_last else al


def get_univariate_df(vartype, feature_data):
    def proc_scores(scores):
        if vartype == "continuous":
            return array_to_list(scores, True)
        elif vartype == "categorical":
            return array_to_list(scores, False)
        else:
            assert False, "wrong vartype, dummy: {}".format(vartype)

    return {
        "value": feature_data["names"],
        "score": proc_scores(feature_data["scores"]),
    }


def build_scores_array(arr, left_type, right_type):
    # any continuous dimension is going to be one less than the names length,
    # so copy the last row/column to the end in these cases
    if left_type == "continuous":
        arr = np.vstack((arr, arr[-1, :]))
    if right_type == "continuous":
        arr = np.column_stack((arr, arr[:, -1]))
    return arr


def get_bivariate_df(all_features, ebm_global, name, data_index):
    def lookup_feature_type(name):
        return all_features[name][0]

    feature_data = ebm_global.data(data_index)
    # left is first dimension, right is second
    left_name, right_name = tuple(name.split(" x "))

    left_type = lookup_feature_type(left_name)
    right_type = lookup_feature_type(right_name)

    left_index = pd.Index(feature_data["left_names"], name="left_value")
    right_index = pd.Index(feature_data["right_names"], name="right_value")

    arr = build_scores_array(feature_data["scores"], left_type, right_type)
    df = (
        pd.DataFrame(arr, index=left_index, columns=right_index)
        .stack()
        .rename("score")
        .reset_index()
        .to_dict(orient="list")
    )

    return {
        "left": {"name": left_name, "type": left_type},
        "right": {"name": right_name, "type": right_type},
        "df": df,
    }


def get_global_scores(ebm_global):
    glob = ebm_global.data()
    return {"variable": glob["names"], "score": glob["scores"]}


def get_univariate_list(ebm_global, all_features):
    return [
        {
            "name": name,
            "vartype": vartype,
            "df": get_univariate_df(vartype, ebm_global.data(i)),
        }
        for name, (vartype, i) in all_features.items()
        if vartype in ["continuous", "categorical"]
    ]


def get_bivariate_list(ebm_global, all_features):
    return [
        get_bivariate_df(all_features, ebm_global, name, i)
        for name, (vartype, i) in all_features.items()
        if vartype == "interaction"
    ]


def get_model_dict(ebm):
    ebm_global = ebm.explain_global()
    all_features = {
        n: (t, i)
        for i, (n, t) in enumerate(
            map(tuple, ebm_global.selector[["Name", "Type"]].to_numpy())
        )
    }
    return {
        "global": get_global_scores(ebm_global),
        "univariate": get_univariate_list(ebm_global, all_features),
        "bivariate": get_bivariate_list(ebm_global, all_features),
    }


def write_predictions(ebm, X_test, y_test):
    y_pred = pd.DataFrame(
        {
            "prob": ebm.predict_proba(X_test)[::, 1],
            "label": y_test["label"],
        }
    )
    y_pred.to_csv(snakemake.output["predictions"], index=False)


def write_model_json(ebm):
    with open(snakemake.output["model"], "w") as f:
        json.dump(get_model_dict(ebm), f)


def main():
    ebm = read_pickle(snakemake.input["model"])
    X_test = pd.read_csv(snakemake.input["test_x"])
    y_test = pd.read_csv(snakemake.input["test_y"])
    write_predictions(ebm, X_test, y_test)
    write_model_json(ebm)


main()
