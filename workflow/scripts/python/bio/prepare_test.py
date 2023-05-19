import pandas as pd
from typing import Any
from common.tsv import write_tsv
from common.io import setup_logging
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
    filter_col = sconf.feature_definitions.vcf.filter
    label_col = sconf.feature_definitions.label_name
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        [cfg.FeatureKey(c) for c in cfg.IDX_COLS],
        filter_col,
        cfg.FeatureKey(label_col),
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
        [cfg.FeatureKey(c) for c in cfg.IDX_COLS],
        df,
    )
    write_tsv(xpath, processed)


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    sin = smk.input
    sout = smk.output
    wcs = smk.wildcards
    variables = sconf.testkey_to_variables(
        cfg.ModelKey(wcs["model_key"]),
        cfg.TestKey(wcs["test_key"]),
    )
    df = pd.read_table(sin["annotated"][0]).assign(
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
