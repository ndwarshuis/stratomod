import json
import pandas as pd
from typing import Any
from os.path import dirname, basename
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
import common.config as cfg
from common.prepare import process_labeled_data, process_unlabeled_data

logger = setup_logging(snakemake.log[0])  # type: ignore


# TODO get rid of this since I don't use vcf input anymore
def read_vcf_input(path: str, input_key: str) -> pd.DataFrame:
    with open(path, "r") as f:
        mapping = {basename(dirname(k)): v for k, v in json.load(f).items()}
        return mapping[input_key]


def read_input(
    df_path: str,
    mapping_path: str,
    input_key: str,
    input_col: str,
) -> pd.DataFrame:
    vcf_input = read_vcf_input(mapping_path, input_key)
    return read_tsv(df_path).assign(**{input_col: vcf_input})


def write_labeled(
    xpath: str,
    ypath: str,
    sconf: cfg.StratoMod,
    rconf: cfg.Model,
    filter_col: cfg.FeatureKey,
    df: pd.DataFrame,
) -> None:
    label_col = sconf.feature_meta.label_name
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        sconf.feature_meta.all_index_cols(),
        filter_col,
        label_col,
        df,
    )
    write_tsv(xpath, processed.drop([label_col], axis=1))
    write_tsv(ypath, processed[label_col].to_frame())


def write_unlabeled(
    xpath: str,
    sconf: cfg.StratoMod,
    rconf: cfg.Model,
    df: pd.DataFrame,
) -> None:
    processed = process_unlabeled_data(
        rconf.features,
        sconf.feature_meta.all_index_cols(),
        df,
    )
    write_tsv(xpath, processed)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    sin = smk.input
    sout = smk.output
    wcs = smk.wildcards
    raw_df = read_input(
        sin["annotated"][0],
        sin["paths"],
        wcs["input_key"],
        sconf.feature_meta.vcf.input_name,
    )
    rconf = sconf.models[cfg.ModelKey(wcs.model_key)]
    if "test_y" in dict(sout):
        write_labeled(
            sout["test_x"],
            sout["test_y"],
            sconf,
            rconf,
            sconf.feature_meta.vcf.filter_name,
            raw_df,
        )
    else:
        write_unlabeled(sout["test_x"], sconf, rconf, raw_df)


main(snakemake, snakemake.config)  # type: ignore
