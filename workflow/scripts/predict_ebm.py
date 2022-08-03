import pandas as pd
from common.cli import setup_logging
from common.ebm import read_model

setup_logging(snakemake.log[0])


def write_csv(path, df):
    df.to_csv(path, header=True, index=False)


def predict_from_x(ebm, df):
    probs, explanations = ebm.predict_and_contrib(df)
    # TODO add column names so a human to read these; likely this will entail
    # getting a list of features as done for postprocessing and training but
    # will also need to add interactions as necessary
    return pd.DataFrame(probs), pd.DataFrame(explanations)


def main():
    sin = snakemake.input
    sout = snakemake.output
    ebm = read_model(sin["model"])
    predict_x = pd.read_csv(sin["predict_x"])
    ps, xs = predict_from_x(ebm, predict_x)
    write_csv(sout["train_x"], ps)
    write_csv(sout["train_y"], xs)


main()
