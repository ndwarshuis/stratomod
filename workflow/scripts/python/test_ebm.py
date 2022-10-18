import pandas as pd
from common.cli import setup_logging
from common.tsv import read_tsv, write_tsv
from common.ebm import read_model
from common.config import lookup_bed_cols_ordered

setup_logging(snakemake.log[0])


def _write_tsv(path, df):
    write_tsv(path, df, header=True)


def predict_from_x(ebm, df):
    probs, explanations = ebm.predict_and_contrib(df)
    return pd.DataFrame(probs), pd.DataFrame(explanations, columns=ebm.feature_names)


def main():
    sin = snakemake.input
    sout = snakemake.output
    ebm = read_model(sin["model"])
    predict_x = read_tsv(sin["test_x"]).drop(
        columns=lookup_bed_cols_ordered(snakemake.config)
    )
    ps, xs = predict_from_x(ebm, predict_x)
    _write_tsv(sout["predictions"], ps)
    _write_tsv(sout["explanations"], xs)


main()
