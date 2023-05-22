import pandas as pd
from typing import Any
from common.io import setup_logging
from common.tsv import write_tsv
from common.ebm import read_model
import common.config as cfg
from interpret.glassbox import ExplainableBoostingClassifier  # type: ignore

setup_logging(snakemake.log[0])  # type: ignore


def _write_tsv(path: str, df: pd.DataFrame) -> None:
    write_tsv(path, df, header=True)


def predict_from_x(
    ebm: ExplainableBoostingClassifier,
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    probs, explanations = ebm.predict_and_contrib(df)
    return pd.DataFrame(probs), pd.DataFrame(explanations, columns=ebm.feature_names)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    sin = smk.input
    sout = smk.output
    ebm = read_model(sin["model"])
    predict_x = pd.read_table(sin["test_x"]).drop(columns=cfg.IDX_COLS)
    ps, xs = predict_from_x(ebm, predict_x)
    _write_tsv(sout["predictions"], ps)
    _write_tsv(sout["explanations"], xs)


main(snakemake, snakemake.config)  # type: ignore
