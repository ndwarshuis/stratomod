import json
import pandas as pd
from functools import partial
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.config import fmt_vcf_feature, lookup_ebm_run
from common.prepare import process_labeled_data

logger = setup_logging(snakemake.log[0])


def read_inputs(paths, input_col):
    eps = [*enumerate(paths)]
    return (
        pd.concat([read_tsv(p).assign(**{input_col: i}) for i, p in eps]),
        # TODO juse use the parent basename here? This should be enough to
        # identify the file (assuming that it is clear which are for indels and
        # snps)
        {p: i for i, p in eps},
    )


def main():
    sconf = snakemake.config
    sout = snakemake.output
    rconf = lookup_ebm_run(sconf, snakemake.wildcards.run_key)
    fconf = sconf["features"]
    label_col = fconf["label"]
    chrom_col = fconf["index"]["chr"]
    _fmt_vcf_feature = partial(fmt_vcf_feature, sconf)
    filter_col = _fmt_vcf_feature("filter")
    input_col = _fmt_vcf_feature("input")
    raw_df, mapped_paths = read_inputs(snakemake.input, input_col)
    processed = process_labeled_data(
        rconf["features"],
        rconf["error_labels"],
        rconf["filtered_are_candidates"],
        chrom_col,
        filter_col,
        label_col,
        raw_df,
    )
    with open(sout["paths"], "w") as f:
        json.dump(mapped_paths, f)
    write_tsv(sout["df"], processed)


main()
