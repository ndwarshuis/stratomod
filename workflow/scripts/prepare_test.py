import json
from os.path import dirname, basename
from functools import partial
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
from common.config import fmt_vcf_feature, lookup_ebm_run
from common.prepare import process_data

logger = setup_logging(snakemake.log[0])


def read_vcf_input(path, input_key):
    with open(path, "r") as f:
        mapping = {basename(dirname(k)): v for k, v in json.load(f).items()}
        return mapping[input_key]


def read_input(df_path, mapping_path, input_key, input_col):
    vcf_input = read_vcf_input(mapping_path, input_key)
    return read_tsv(df_path).assign(**{input_col: vcf_input})


def main():
    sin = snakemake.input
    sout = snakemake.output
    sconf = snakemake.config
    wcs = snakemake.wildcards
    _fmt_vcf_feature = partial(fmt_vcf_feature, sconf)
    raw_df = read_input(
        sin["annotated"],
        sin["paths"],
        wcs["input_key"],
        _fmt_vcf_feature("input"),
    )
    rconf = lookup_ebm_run(sconf, wcs.run_key)
    fconf = sconf["features"]
    label_col = fconf["label"]
    processed = process_data(
        rconf["features"],
        rconf["error_labels"],
        rconf["filtered_are_candidates"],
        fconf["index"]["chr"],
        _fmt_vcf_feature("filter"),
        label_col,
        raw_df,
    )
    write_tsv(sout["test_x"], processed.drop([label_col]))
    write_tsv(sout["test_y"], processed[label_col])


main()
