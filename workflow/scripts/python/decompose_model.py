import json
import pandas as pd
import numpy as np
from common.cli import setup_logging
from common.ebm import read_model

setup_logging(snakemake.log[0])


# TODO there is no reason this can't be done immediately after training
# just to avoid the pickle thing


def array_to_list(arr, repeat_last):
    al = arr.tolist()
    return al + [al[-1]] if repeat_last else al


def get_univariate_df(vartype, feature_data, stdev):
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
        # For some reason, the standard deviations array has an extra 0 in
        # in the front and thus is one longer than the scores array.
        "stdev": proc_scores(stdev[1:]),
    }


def build_scores_array(arr, left_type, right_type):
    # any continuous dimension is going to be one less than the names length,
    # so copy the last row/column to the end in these cases
    if left_type == "continuous":
        arr = np.vstack((arr, arr[-1, :]))
    if right_type == "continuous":
        arr = np.column_stack((arr, arr[:, -1]))
    return arr


def get_bivariate_df(all_features, ebm_global, name, data_index, stdevs):
    def lookup_feature_type(name):
        return all_features[name][0]

    feature_data = ebm_global.data(data_index)
    # left is first dimension, right is second
    left_name, right_name = tuple(name.split(" x "))

    left_type = lookup_feature_type(left_name)
    right_type = lookup_feature_type(right_name)

    left_index = pd.Index(feature_data["left_names"], name="left_value")
    right_index = pd.Index(feature_data["right_names"], name="right_value")

    def stack_array(arr, name):
        return (
            pd.DataFrame(
                build_scores_array(arr, left_type, right_type),
                index=left_index,
                columns=right_index,
            )
            .stack()
            .rename(name)
        )

    # the standard deviations are in an array that has 1 larger shape than the
    # scores array in both directions where the first row/column is all zeros.
    # Not sure why it is all zeros, but in order to make it line up with the
    # scores array we need to shave off the first row/column.
    return {
        "left": {"name": left_name, "type": left_type},
        "right": {"name": right_name, "type": right_type},
        "df": pd.concat(
            [
                stack_array(feature_data["scores"], "score"),
                stack_array(stdevs[data_index][1:, 1:], "stdev"),
            ],
            axis=1,
        )
        .reset_index()
        .to_dict(orient="list"),
    }


def get_global_scores(ebm_global):
    glob = ebm_global.data()
    return {"variable": glob["names"], "score": glob["scores"]}


def get_univariate_list(ebm_global, all_features, stdevs):
    return [
        {
            "name": name,
            "vartype": vartype,
            "df": get_univariate_df(vartype, ebm_global.data(i), stdevs[i]),
        }
        for name, (vartype, i) in all_features.items()
        if vartype in ["continuous", "categorical"]
    ]


def get_bivariate_list(ebm_global, all_features, stdevs):
    return [
        get_bivariate_df(all_features, ebm_global, name, i, stdevs)
        for name, (vartype, i) in all_features.items()
        if vartype == "interaction"
    ]


def get_model_dict(ebm):
    ebm_global = ebm.explain_global()
    stdevs = ebm.term_standard_deviations_
    all_features = {
        n: (t, i)
        for i, (n, t) in enumerate(
            map(tuple, ebm_global.selector[["Name", "Type"]].to_numpy())
        )
    }
    return {
        "global": get_global_scores(ebm_global),
        "intercept": ebm.intercept_[0],
        "univariate": get_univariate_list(ebm_global, all_features, stdevs),
        "bivariate": get_bivariate_list(ebm_global, all_features, stdevs),
    }


def write_predictions(ebm, X_test, y_test, label):
    y_pred = pd.DataFrame(
        {
            "prob": ebm.predict_proba(X_test)[::, 1],
            "label": y_test[label],
        }
    )
    y_pred.to_csv(snakemake.output["predictions"], index=False)


def write_model_json(ebm):
    with open(snakemake.output["model"], "w") as f:
        json.dump(get_model_dict(ebm), f)


def main():
    ebm = read_model(snakemake.input["model"])
    label = snakemake.config["features"]["label"]
    X_test = pd.read_csv(snakemake.input["test_x"])
    y_test = pd.read_csv(snakemake.input["test_y"])
    write_predictions(ebm, X_test, y_test, label)
    write_model_json(ebm)


main()
