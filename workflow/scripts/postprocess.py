import json
import pandas as pd
import numpy as np
from more_itertools import duplicates_everseen
from common.tsv import read_tsv, write_tsv
from common.cli import make_io_parser, add_config_arg


def make_parser():
    parser = make_io_parser(
        "transform and impute dataframe",
        "the input dataframe",
        "the processed dataframe",
    )
    add_config_arg("the configuration to use for processing", parser)
    return parser


def process_series(opts, ser):
    log_trans = opts["log_transform"]
    fillval = opts["fill_na"]
    _ser = pd.to_numeric(ser, errors="coerce")
    return (np.log(_ser) if log_trans else _ser).fillna(fillval)


def check_columns(wanted_cols, df_cols):
    wanted_set = set(wanted_cols)
    df_set = set(df_cols)
    assert len(wanted_set) == len(
        wanted_cols
    ), "duplicate configuration features detected: {}".format(
        list(duplicates_everseen(wanted_cols)),
    )
    assert len(df_set) == len(
        df_cols
    ), "input dataframe has duplicate columns: {}".format(
        list(duplicates_everseen(df_cols)),
    )
    assert (
        df_set >= wanted_set
    ), "configuration features must be a subset of columns in input dataframe"


def select_columns(config, df):
    # TODO don't hardcode this
    label = "label"
    wanted_cols = list(config) + [label]
    check_columns(wanted_cols, df.columns.tolist())
    return df[wanted_cols]


def process_data(config, df):
    for col, opts in config.items():
        df[col] = process_series(opts, df[col])
    # select columns after transforms to avoid pandas asking me to make a
    # deep copy (which will happen on a slice of a slice)
    return select_columns(config, df)


def main():
    args = make_parser().parse_args()
    config = json.loads(args.config)
    raw = read_tsv(args.input)
    processed = process_data(config, raw)
    write_tsv(args.output, processed)


main()
