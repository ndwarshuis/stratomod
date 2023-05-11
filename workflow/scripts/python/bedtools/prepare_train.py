import pandas as pd
import common.config as cfg
from typing import Any
from common.tsv import write_tsv
from common.cli import setup_logging
from common.prepare import process_labeled_data

logger = setup_logging(snakemake.log[0])  # type: ignore


def read_query(
    config: cfg.StratoMod, path: str, key: cfg.LabeledQueryKey
) -> pd.DataFrame:
    variables = config.querykey_to_variables(key)
    return pd.read_table(path).assign(**{str(k): v for k, v in variables.items()})


def read_queries(
    config: cfg.StratoMod,
    paths: dict[cfg.LabeledQueryKey, str],
) -> pd.DataFrame:
    # TODO this is weird, why do I need the [0] here?
    return pd.concat([read_query(config, path[0], key) for key, path in paths.items()])


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    rconf = sconf.models[cfg.ModelKey(cfg.ModelKey(smk.wildcards.model_key))]
    fconf = sconf.feature_names
    raw_df = read_queries(sconf, smk.input)
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        [cfg.FeatureKey(x) for x in fconf.all_index_cols()],
        cfg.FeatureKey(fconf.vcf.fmt_name(lambda x: x.filter)),
        cfg.FeatureKey(fconf.label_name),
        raw_df,
    )
    write_tsv(smk.output["df"], processed)


main(snakemake, snakemake.config)  # type: ignore
