import pandas as pd
import common.config as cfg
from typing import Dict, Any
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.prepare import process_labeled_data

logger = setup_logging(snakemake.log[0])  # type: ignore


def read_query(
    config: cfg.StratoMod, path: str, key: cfg.LabeledQueryKey
) -> pd.DataFrame:
    variables = config.labeled_queries[key].variables
    return read_tsv(path).assign(**{str(k): v for k, v in variables.items()})


def read_queries(
    config: cfg.StratoMod,
    paths: Dict[cfg.LabeledQueryKey, str],
) -> pd.DataFrame:
    return pd.concat([read_query(config, path, key) for key, path in paths.items()])


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    rconf = sconf.models[cfg.ModelKey(cfg.ModelKey(smk.wildcards.model_key))]
    fconf = sconf.feature_meta
    raw_df = read_queries(sconf, smk.input)
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        sconf.feature_meta.all_index_cols(),
        fconf.vcf.filter_name,
        fconf.label_name,
        raw_df,
    )
    write_tsv(smk.output["df"], processed)


main(snakemake, snakemake.config)  # type: ignore
