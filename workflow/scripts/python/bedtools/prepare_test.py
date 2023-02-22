import pandas as pd
from typing import Any
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
import common.config as cfg
from common.prepare import process_labeled_data, process_unlabeled_data

logger = setup_logging(snakemake.log[0])  # type: ignore


def write_labeled(
    xpath: str,
    ypath: str,
    sconf: cfg.StratoMod,
    rconf: cfg.Model,
    df: pd.DataFrame,
) -> None:
    filter_col = sconf.feature_names.vcf.filter_name
    label_col = sconf.feature_names.label_name
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        sconf.feature_names.all_index_cols(),
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
        sconf.feature_names.all_index_cols(),
        df,
    )
    write_tsv(xpath, processed)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    sin = smk.input
    sout = smk.output
    wcs = smk.wildcards
    variables = sconf.testkey_to_variables(
        cfg.ModelKey(wcs["model_key"]),
        cfg.RunKey(wcs["run_key"]),
        cfg.TestKey(wcs["test_key"]),
    )
    df = read_tsv(sin["annotated"][0]).assign(
        **{str(k): v for k, v in variables.items()}
    )
    rconf = sconf.models[cfg.ModelKey(wcs.model_key)]
    if "test_y" in dict(sout):
        write_labeled(
            sout["test_x"],
            sout["test_y"],
            sconf,
            rconf,
            df,
        )
    else:
        write_unlabeled(sout["test_x"], sconf, rconf, df)


main(snakemake, snakemake.config)  # type: ignore
