import pandas as pd
import yaml
from typing import Any
from more_itertools import flatten
from sklearn.model_selection import train_test_split  # type: ignore
from interpret.glassbox import ExplainableBoostingClassifier  # type: ignore
from common.tsv import write_tsv
from common.cli import setup_logging
from common.ebm import write_model
import common.config as cfg

logger = setup_logging(snakemake.log[0])  # type: ignore


def _write_tsv(smk: Any, key: str, df: pd.DataFrame) -> None:
    write_tsv(smk.output[key], df, header=True)


def dump_config(smk: Any, config: cfg.Model) -> None:
    with open(smk.output["config"], "w") as f:
        yaml.dump(config, f)


def get_interactions(
    df_columns: list[cfg.FeatureKey],
    iconfig: int | cfg.InteractionSpec,
) -> int | list[list[int]]:
    def expand_interactions(i: cfg.InteractionSpec_) -> list[list[int]]:
        if isinstance(i, str):
            return [
                [df_columns.index(i), c] for c, f in enumerate(df_columns) if f != i
            ]
        else:
            return [[df_columns.index(i.f1), df_columns.index(i.f2)]]

    if isinstance(iconfig, int):
        return iconfig
    else:
        return [*flatten(expand_interactions(i) for i in iconfig)]


def train_ebm(
    smk: Any,
    sconf: cfg.StratoMod,
    rconf: cfg.Model,
    df: pd.DataFrame,
) -> None:
    label = sconf.feature_names.label

    def strip_coords(df: pd.DataFrame) -> pd.DataFrame:
        return df.drop(columns=sconf.feature_names.all_index_cols())

    features = rconf.features
    feature_names = [
        k if v.alt_name is None else v.alt_name for k, v in features.items()
    ]
    misc_params = rconf.ebm_settings.misc_parameters

    if misc_params.downsample is not None:
        df = df.sample(frac=misc_params.downsample)

    train_cols = [c for c in df.columns if c != label]
    X = df[train_cols]
    y = df[label]

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        **rconf.ebm_settings.split_parameters,
    )

    ebm_config = rconf.ebm_settings.classifier_parameters

    cores = smk.threads

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
        feature_types=[f.feature_type for f in features.values()],
        interactions=get_interactions(feature_names, rconf.interactions),
        n_jobs=cores,
        **ebm_config,
    )
    ebm.fit(strip_coords(X_train), y_train)

    write_model(smk.output["model"], ebm)
    _write_tsv(smk, "train_x", X_train)
    _write_tsv(smk, "train_y", y_train)
    _write_tsv(smk, "test_x", X_test)
    _write_tsv(smk, "test_y", y_test)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    rconf = sconf.models[smk.wildcards.model_key]
    df = pd.read_table(smk.input[0])
    train_ebm(smk, sconf, rconf, df)
    dump_config(smk, rconf)


main(snakemake, snakemake.config)  # type: ignore
