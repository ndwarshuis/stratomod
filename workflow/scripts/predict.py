# import pandas as pd
# import numpy as np
from common.cli import setup_logging
from common.ebm import read_model

setup_logging(snakemake.log[0])


def main():
    ebm = read_model(snakemake.input["model"])
    predict_X = pd.read_csv(snakemake.input["predict_x"])
    # TODO write these (obviously)

    # TODO should add feature names to the top of the explanations so we can
    # actually interpret them
    probs, explanations = ebm.predict_and_contrib(predict_X)


main()
