import argparse
import pickle
import json
import yaml
from os.path import join
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider
import pandas as pd
from sklearn.model_selection import train_test_split
from interpret.glassbox import ExplainableBoostingClassifier
from common.tsv import read_tsv
from common.cli import add_input_arg, add_config_arg

MODEL_FILE = "model.pickle"
TRAIN_X_FILE = "train_x.pickle"
TEST_X_FILE = "test_x.pickle"
TRAIN_Y_FILE = "train_y.pickle"
TEST_Y_FILE = "test_y.pickle"
CONFIG_FILE = "config.yml"


def make_parser():
    parser = argparse.ArgumentParser(description="Train an EBM model.")
    add_input_arg("the training dataframe", parser)
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="the output directory in which the trained results will be written",
    )
    add_config_arg("the configuration to use for training", parser)
    return parser


def write_pickle(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def train_ebm(config, label, outdir, df):
    ebm_config = config["settings"]
    # TODO don't hardcode this
    # df = df.sample(frac=0.01)
    train_cols = [c for c in df.columns if c != label]
    X = df[train_cols]
    y = df[label]

    # TODO don't hardcode this
    seed = 1

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=ebm_config["test_size"],
        random_state=seed,
    )

    ebm = ExplainableBoostingClassifier(random_state=seed)
    ebm.fit(X_train, y_train)

    write_pickle(join(outdir, MODEL_FILE), ebm)
    write_pickle(join(outdir, TRAIN_X_FILE), X_train)
    write_pickle(join(outdir, TRAIN_Y_FILE), y_train)
    write_pickle(join(outdir, TEST_X_FILE), X_test)
    write_pickle(join(outdir, TEST_Y_FILE), y_test)


def dump_config(config, outdir):
    with open(join(outdir, CONFIG_FILE), "w") as f:
        yaml.dump(config, f)


def main():
    args = make_parser().parse_args()
    config = json.loads(args.config)
    df = read_tsv(args.input)
    # TODO don't hardcode the label
    train_ebm(config, "label", args.output_dir, df)
    dump_config(config, args.output_dir)


main()
