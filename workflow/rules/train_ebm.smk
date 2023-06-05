import re
import json
import sys
import os
from pathlib import Path
import subprocess as sp
from collections import namedtuple
from more_itertools import flatten, partition, unzip
import scripts.python.common.config as cfg

# giant hack to allow sourcing other files in rmd scripts (seems like this
# should be built-in...oh well)

_python_path = next(
    (
        p
        for p in map(lambda p: Path(os.getcwd()) / p, sys.path)
        if p.is_absolute() and p.exists() and p.name == "python"
    ),
    None,
)
assert _python_path is not None
rmd_lib_path = _python_path.parent / "rmarkdown" / "common"


################################################################################
# annotate variants with features


annotated_tsv = cfg.wildcard_ext("vartype_key", "tsv.gz")
annotated_log = cfg.wildcard_ext("vartype_key", "log")
annotated_bench = cfg.wildcard_ext("vartype_key", "bench")


def annotation_input(tsv_path, query_key):
    rsk = config.querykey_to_refsetkey(query_key)
    rk = config.refsetkey_to_refkey(rsk)
    return {
        "variants": expand(tsv_path, allow_missing=True, refset_key=rsk),
        "features": [
            *expand(rmsk_targets(rk), refset_key=rsk),
            *expand(
                [
                    *rules.get_tandem_repeats.output,
                    rules.get_mappability.output.high,
                    rules.get_mappability.output.low,
                    *rules.get_segdups.output,
                ],
                allow_missing=True,
                refset_key=rsk,
            ),
            *expand(
                rules.get_homopolymers.output,
                allow_missing=True,
                base=cfg.Base.all(),
                refset_key=rsk,
            ),
        ],
    }


rule annotate_labeled_variants:
    input:
        unpack(
            lambda w: annotation_input(rules.concat_labeled_tsvs.output, w.l_query_key)
        ),
    output:
        config.annotated_res_dir(labeled=True, log=False) / annotated_tsv,
    conda:
        "../envs/bio.yml"
    log:
        config.annotated_res_dir(labeled=True, log=True) / annotated_log,
    benchmark:
        config.annotated_res_dir(labeled=True, log=True) / annotated_bench
    script:
        "../scripts/python/bio/annotate_variants.py"


use rule annotate_labeled_variants as annotate_unlabeled_variants with:
    input:
        unpack(
            lambda wildcards: annotation_input(
                rules.parse_unlabeled_vcf.output, wildcards.ul_query_key
            )
        ),
    output:
        config.annotated_res_dir(labeled=False, log=False) / annotated_tsv,
    log:
        config.annotated_res_dir(labeled=False, log=True) / annotated_log,
    benchmark:
        config.annotated_res_dir(labeled=False, log=True) / annotated_bench


################################################################################
# summarize annotated input

summary_output = cfg.wildcard_format("{}_summary.html", "vartype_key")
summary_bench = cfg.wildcard_format("{}_summary.bench", "vartype_key")


rule summarize_labeled_annotated_variants:
    input:
        rules.annotate_labeled_variants.output,
    output:
        config.annotated_res_dir(labeled=True, log=False) / summary_output,
    conda:
        "../envs/summary.yml"
    benchmark:
        config.annotated_res_dir(labeled=True, log=True) / summary_bench
    params:
        label_col=config.feature_definitions.label,
        columns=config.feature_definitions.non_summary_cols,
        query_key=lambda w: w.l_query_key,
        lib_path=rmd_lib_path,
    script:
        "../scripts/rmarkdown/summary/input_summary.Rmd"


use rule summarize_labeled_annotated_variants as summarize_unlabeled_annotated_variants with:
    input:
        rules.annotate_unlabeled_variants.output,
    output:
        config.annotated_res_dir(labeled=False, log=False) / summary_output,
    benchmark:
        config.annotated_res_dir(labeled=False, log=True) / summary_bench
    params:
        label_col=None,
        columns=config.feature_definitions.non_summary_cols,
        query_key=lambda w: w.ul_query_key,


################################################################################
# training
#
# assume that this will take care of test/train split, actual training, and
# pickling


rule prepare_train_data:
    input:
        # NOTE: name these inputs according to their labeled_query_keys so
        # they can be used in the script along with each file path
        unpack(
            lambda w: {
                k: expand(
                    rules.annotate_labeled_variants.output,
                    allow_missing=True,
                    l_query_key=k,
                )
                for k in config.modelkey_to_train_querykeys(w.model_key)
            }
        ),
    output:
        df=config.model_train_res_dir(log=False) / "train.tsv.gz",
    log:
        config.model_train_res_dir(log=True) / "prepare.log",
    benchmark:
        config.model_train_res_dir(log=True) / "prepare.bench"
    conda:
        "../envs/bio.yml"
    script:
        "../scripts/python/bio/prepare_train.py"


rule train_model:
    input:
        rules.prepare_train_data.output,
    output:
        **{
            n: str((config.model_train_res_dir(log=False) / n).with_suffix(".tsv.gz"))
            for n in [
                "train_x",
                "train_y",
                "test_x",
                "test_y",
            ]
        },
        model=config.model_train_res_dir(log=False) / "model.pickle",
        config=config.model_train_res_dir(log=False) / "config.yml",
    conda:
        "../envs/ebm.yml"
    log:
        config.model_train_res_dir(log=True) / "model.log",
    threads: 1
    benchmark:
        config.model_train_res_dir(log=True) / "model.bench"
    script:
        "../scripts/python/ebm/train_ebm.py"


rule decompose_model:
    input:
        **rules.train_model.output,
    output:
        model=config.model_train_res_dir(log=False) / "model.json",
        predictions=config.model_train_res_dir(log=False) / "predictions.tsv.gz",
        train_predictions=config.model_train_res_dir(log=False)
        / "train_predictions.tsv.gz",
    conda:
        "../envs/ebm.yml"
    log:
        config.model_train_res_dir(log=True) / "decompose.log",
    script:
        "../scripts/python/ebm/decompose_model.py"


rule summarize_model:
    input:
        **rules.decompose_model.output,
        train_x=rules.train_model.output.train_x,
        train_y=rules.train_model.output.train_y,
    output:
        config.model_train_res_dir(log=False) / "summary.html",
    conda:
        "../envs/summary.yml"
    benchmark:
        config.model_train_res_dir(log=True) / "summary.bench"
    params:
        # TODO this is because we can't convert enums yet
        features=lambda w: {
            k: v.r_dict for k, v in config.models[w.model_key].features.items()
        },
        error_labels=lambda w: [
            x.value for x in config.models[w.model_key].error_labels
        ],
        lib_path=rmd_lib_path,
    script:
        "../scripts/rmarkdown/summary/train_summary.Rmd"


################################################################################
# prepare test data

prepare_x_file = "test_x.tsv.gz"
prepare_log_file = "prepare.log"
prepare_bench_file = "prepare.bench"


def test_data_input(annotated_path, key, wildcards):
    return {
        "annotated": expand(
            annotated_path,
            vartype_key=wildcards.vartype_key,
            **{
                key: config.testkey_to_querykey(
                    wildcards.model_key,
                    wildcards.test_key,
                )
            },
        ),
    }


# NOTE: the script will treat the input differently depending on if we want
# test_y (in this case it will perform all the label-related logic)
rule prepare_labeled_test_data:
    input:
        unpack(
            partial(
                test_data_input,
                rules.annotate_labeled_variants.output,
                "l_query_key",
            )
        ),
    output:
        test_x=config.model_test_res_dir(labeled=True, log=False) / prepare_x_file,
        test_y=config.model_test_res_dir(labeled=True, log=False) / "test_y.tsv.gz",
    log:
        config.model_test_res_dir(labeled=True, log=True) / prepare_log_file,
    benchmark:
        config.model_test_res_dir(labeled=True, log=True) / prepare_bench_file
    conda:
        "../envs/bio.yml"
    script:
        "../scripts/python/bio/prepare_test.py"


use rule prepare_labeled_test_data as prepare_unlabeled_test_data with:
    input:
        unpack(
            partial(
                test_data_input,
                rules.annotate_unlabeled_variants.output,
                "ul_query_key",
            )
        ),
    output:
        test_x=config.model_test_res_dir(labeled=False, log=False) / prepare_x_file,
    log:
        config.model_test_res_dir(labeled=False, log=True) / prepare_log_file,
    benchmark:
        config.model_test_res_dir(labeled=False, log=True) / prepare_bench_file


################################################################################
# test model


test_bench_file = "test.bench"
test_log_file = "test.log"


def test_ebm_input(x_path):
    return {"model": rules.train_model.output.model, "test_x": x_path}


def test_ebm_output(test_path):
    return {
        "predictions": test_path / "predictions.tsv.gz",
        "explanations": test_path / "explanations.tsv.gz",
    }


rule test_labeled_ebm:
    input:
        **test_ebm_input(rules.prepare_labeled_test_data.output.test_x),
    output:
        **test_ebm_output(config.model_test_res_dir(labeled=True, log=False)),
    conda:
        "../envs/ebm.yml"
    log:
        config.model_test_res_dir(labeled=True, log=True) / test_log_file,
    benchmark:
        config.model_test_res_dir(labeled=True, log=True) / test_bench_file
    script:
        "../scripts/python/ebm/test_ebm.py"


use rule test_labeled_ebm as test_unlabeled_ebm with:
    input:
        **test_ebm_input(rules.prepare_unlabeled_test_data.output.test_x),
    output:
        **test_ebm_output(config.model_test_res_dir(labeled=False, log=False)),
    log:
        config.model_test_res_dir(labeled=False, log=True) / test_log_file,
    benchmark:
        config.model_test_res_dir(labeled=False, log=True) / test_bench_file


################################################################################
# summarize test

test_summary_file = "summary.html"
test_summary_bench = "summary.bench"


rule summarize_labeled_test:
    input:
        **rules.test_labeled_ebm.output,
        truth_y=rules.prepare_labeled_test_data.output.test_y,
    output:
        config.model_test_res_dir(labeled=True, log=False) / test_summary_file,
    conda:
        "../envs/summary.yml"
    benchmark:
        config.model_test_res_dir(labeled=True, log=True) / test_summary_bench
    params:
        query_key=lambda w: config.testkey_to_querykey(w.model_key, w.test_key),
        lib_path=rmd_lib_path,
    script:
        "../scripts/rmarkdown/summary/test_summary.Rmd"


use rule summarize_labeled_test as summarize_unlabeled_test with:
    input:
        **rules.test_unlabeled_ebm.output,
    output:
        config.model_test_res_dir(labeled=False, log=False) / test_summary_file,
    benchmark:
        config.model_test_res_dir(labeled=False, log=True) / test_summary_bench


################################################################################
# global targets


rule all_summary:
    input:
        config.summary_targets(
            rules.summarize_labeled_annotated_variants.output,
            rules.summarize_unlabeled_annotated_variants.output,
        ),


rule all_ebm:
    input:
        config.model_targets(
            rules.summarize_model.output,
            rules.summarize_labeled_test.output,
            rules.summarize_unlabeled_test.output,
        ),
