import pandas as pd
from common.cli import setup_logging
from common.tsv import read_tsv
from common.ebm import read_model

setup_logging(snakemake.log[0])


def write_csv(path, df):
    df.to_csv(path, header=True, index=False)


def predict_from_x(ebm, df):
    probs, explanations = ebm.predict_and_contrib(df)
    return pd.DataFrame(probs), pd.DataFrame(explanations, columns=ebm.feature_names)


def main():
    sin = snakemake.input
    sout = snakemake.output
    ebm = read_model(sin["model"])
    predict_x = read_tsv(sin["test_x"])
    ps, xs = predict_from_x(ebm, predict_x)
    write_csv(sout["predictions"], ps)
    write_csv(sout["explanations"], xs)


main()
