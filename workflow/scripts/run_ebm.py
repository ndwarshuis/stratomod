import random
import pickle
import yaml
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider
from sklearn.model_selection import train_test_split
from interpret.glassbox import ExplainableBoostingClassifier
from common.tsv import read_tsv


def write_model(obj):
    with open(snakemake.output["model"], "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def write_csv(key, df):
    df.to_csv(snakemake.output[key], header=True, index=False)


def dump_config(config):
    with open(snakemake.output["config"], "w") as f:
        yaml.dump(config, f)


def get_interactions(df_columns, iconfig):
    # ASSUME type is int | [[int]]
    if isinstance(iconfig, int):
        return iconfig
    return [[df_columns.index(n) for n in names] for names in iconfig]


def train_ebm(config, label, df):
    features = config["features"]
    misc_params = config["ebm_settings"]["misc_parameters"]

    if not misc_params["downsample"] is None:
        df = df.sample(frac=misc_params["downsample"])

    train_cols = [c for c in df.columns if c != label]
    X = df[train_cols]
    y = df[label]

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        **config["ebm_settings"]["split_parameters"],
    )

    ebm_config = config["ebm_settings"]["classifier_parameters"]

    if ebm_config["random_state"] is None:
        ebm_config["random_state"] = random.randrange(0, 420420)

    ebm = ExplainableBoostingClassifier(
        # NOTE the EBM docs show them explicitly adding interactions here like
        # 'F1 x F2' but it appears to work when I specify them separately via
        # the 'interactions' parameter
        feature_names=list(features),
        feature_types=[f["feature_type"] for f in features.values()],
        interactions=get_interactions(list(features), config["interactions"]),
        **ebm_config,
    )
    ebm.fit(X_train, y_train)

    write_model(ebm)
    write_csv("train_x", X_train)
    write_csv("train_y", y_train)
    write_csv("test_x", X_test)
    write_csv("test_y", y_test)


def main():
    params = snakemake.params
    df = read_tsv(snakemake.input[0])
    # TODO don't hardcode the label
    train_ebm(params.config, "label", df)
    dump_config(params.config)


main()
