import json
import pandas as pd
import common.config as cfg
from typing import List, Tuple, Dict
from functools import partial
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.prepare import process_labeled_data

logger = setup_logging(snakemake.log[0])  # type: ignore


def read_inputs(
    paths: List[str],
    input_col: str,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    eps = [*enumerate(paths)]
    return (
        pd.concat([read_tsv(p).assign(**{input_col: i}) for i, p in eps]),
        # TODO just use the parent basename here? This should be enough to
        # identify the file (assuming that it is clear which are for indels and
        # snps)
        {p: i for i, p in eps},
    )


def main(smk, sconf: cfg.StratoMod) -> None:
    sout = smk.output
    rconf = sconf.models[smk.wildcards.model_key]
    fconf = sconf.feature_meta
    label_col = fconf.label
    _fmt_vcf_feature = partial(sconf.feature_meta.vcf.fmt_feature, sconf)
    filter_col = _fmt_vcf_feature("filter")
    input_col = _fmt_vcf_feature("input")
    raw_df, mapped_paths = read_inputs(smk.input, input_col)
    processed = process_labeled_data(
        rconf.features,
        rconf.error_labels,
        rconf.filtered_are_candidates,
        sconf.feature_meta.all_index_cols(),
        filter_col,
        label_col,
        raw_df,
    )
    with open(sout["paths"], "w") as f:
        json.dump(mapped_paths, f)
    write_tsv(sout["df"], processed)


main(snakemake, snakemake.config)  # type: ignore
