import random
import yaml
from more_itertools import flatten
from sklearn.model_selection import train_test_split
from interpret.glassbox import ExplainableBoostingClassifier
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.ebm import write_model
from common.config import lookup_ebm_run

logger = setup_logging(snakemake.log[0])


def _write_tsv(key, df):
    write_tsv(snakemake.output[key], df, header=True)


def dump_config(config):
    with open(snakemake.output["config"], "w") as f:
        yaml.dump(config, f)


def get_interactions(df_columns, iconfig):
    # ASSUME type is int | [str] | [[str]]
    def expand_interactions(i):
        if isinstance(i, str):
            return [
                [df_columns.index(i), c] for c, f in enumerate(df_columns) if f != i
            ]
        else:
            return [[df_columns.index(feature) for feature in i]]

    if isinstance(iconfig, int):
        return iconfig
    else:
        return [*flatten(expand_interactions(i) for i in iconfig)]


def train_ebm(config, label, df):
    features = config["features"]
    feature_names = [
        k if v["alt_name"] is None else v["alt_name"] for k, v in features.items()
    ]
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

    cores = snakemake.threads

    logger.info(
        "Training EBM with %d features and %d cores",
        len(features),
        cores,
    )

    ebm = ExplainableBoostingClassifier(
        # NOTE the EBM docs show them explicitly adding interactions here like
        # 'F1 x F2' but it appears to work when I specify them separately via
        # the 'interactions' parameter
        feature_names=feature_names,
        feature_types=[f["feature_type"] for f in features.values()],
        interactions=get_interactions(feature_names, config["interactions"]),
        n_jobs=cores,
        **ebm_config,
    )
    ebm.fit(X_train, y_train)

    write_model(snakemake.output["model"], ebm)
    _write_tsv("train_x", X_train)
    _write_tsv("train_y", y_train)
    _write_tsv("test_x", X_test)
    _write_tsv("test_y", y_test)


def main():
    sconf = snakemake.config
    rconf = lookup_ebm_run(sconf, snakemake.wildcards.run_key)
    df = read_tsv(snakemake.input[0])
    train_ebm(rconf, sconf["features"]["label"], df)
    dump_config(rconf)


main()
