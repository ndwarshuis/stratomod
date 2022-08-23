import re
import json
import subprocess as sp
from collections import namedtuple
from more_itertools import flatten, partition, unzip
from scripts.python.common.config import attempt_mem_gb

labeled_dir = "labeled"
unlabeled_dir = "unlabeled"


INPUT_DELIM = "&"


################################################################################
# add annotations


annotated_tsv = "{filter_key}.tsv"
annotated_log = "{filter_key}.log"
annotated_bench = "{filter_key}.bench"

annotated_dir = Path("annotated_input")
rel_annotated_labeled_dir = annotated_dir / labeled_dir / "{input_key}"
rel_annotated_unlabeled_dir = annotated_dir / unlabeled_dir / "{input_key}"

labeled_annotated_dir = results_dir / rel_annotated_labeled_dir
unlabeled_annotated_dir = results_dir / rel_annotated_unlabeled_dir


def annotation_input(tsv_path):
    return {
        "variants": tsv_path,
        "annotations": [
            rules.get_repeat_masker_classes.output,
            rules.get_tandem_repeats.output,
            rules.get_mappability.output.high,
            rules.get_mappability.output.low,
            rules.get_segdups.output,
            expand(
                rules.get_homopolymers.output,
                bases=config["features"]["homopolymers"]["bases"],
                allow_missing=True,
            ),
        ],
    }


rule annotate_labeled_tsv:
    input:
        **annotation_input(rules.concat_labeled_tsvs.output),
    output:
        labeled_annotated_dir / annotated_tsv,
    conda:
        envs_path("bedtools.yml")
    log:
        log_dir / rel_annotated_labeled_dir / annotated_log,
    benchmark:
        labeled_annotated_dir / annotated_bench
    resources:
        mem_mb=attempt_mem_gb(32),
    script:
        python_path("annotate.py")


use rule annotate_labeled_tsv as annotate_unlabeled_tsv with:
    input:
        **annotation_input(rules.parse_unlabeled_vcf.output),
    output:
        unlabeled_annotated_dir / annotated_tsv,
    log:
        log_dir / rel_annotated_unlabeled_dir / annotated_log,
    benchmark:
        unlabeled_annotated_dir / annotated_bench


################################################################################
# summarize annotated input

summary_output = "{filter_key}_summary.html"
summary_bench = "{filter_key}_summary.bench"


rule summarize_labeled_input:
    input:
        rules.annotate_labeled_tsv.output,
    output:
        labeled_annotated_dir / summary_output,
    conda:
        envs_path("rmarkdown.yml")
    benchmark:
        labeled_annotated_dir / summary_bench
    resources:
        mem_mb=attempt_mem_gb(16),
    params:
        has_label=True,
    script:
        rmd_path("input_summary.Rmd")


use rule summarize_labeled_input as summarize_unlabeled_input with:
    input:
        rules.annotate_unlabeled_tsv.output,
    output:
        unlabeled_annotated_dir / summary_output,
    benchmark:
        unlabeled_annotated_dir / summary_bench
    params:
        has_label=False,


################################################################################
# training
#
# assume that this will take care of test/train split, actual training, and
# pickling


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


git_tag = get_git_tag()

ebm_dir = Path("ebm") / f"{git_tag}-{{input_keys}}-{{filter_key}}-{{run_key}}"

train_results_dir = results_dir / ebm_dir
train_log_dir = log_dir / ebm_dir


rule prepare_train_data:
    input:
        lambda wildcards: expand(
            rules.annotate_labeled_tsv.output,
            allow_missing=True,
            input_key=wildcards.input_keys.split(INPUT_DELIM),
        ),
    output:
        df=train_results_dir / "train.tsv",
        paths=train_results_dir / "input_paths.json",
    log:
        train_log_dir / "prepare.log",
    benchmark:
        train_log_dir / "prepare.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        python_path("prepare_train.py")


rule train_model:
    input:
        rules.prepare_train_data.output,
    output:
        **{
            n: str((train_results_dir / n).with_suffix(".csv"))
            for n in [
                "train_x",
                "train_y",
                "test_x",
                "test_y",
            ]
        },
        model=train_results_dir / "model.pickle",
        config=train_results_dir / "config.yml",
    conda:
        envs_path("ebm.yml")
    log:
        train_log_dir / "model.log",
    threads: 1
    # ASSUME total memory is proportional the number of EBMs that are trained in
    # parallel (see outer_bags parameter in the EBM function call)
    resources:
        mem_mb=lambda wildcards, threads, attempt: 16000 * threads * 2 ** (attempt - 1),
    benchmark:
        train_results_dir / "model.bench"
    script:
        python_path("train_ebm.py")


rule decompose_model:
    input:
        **rules.train_model.output,
    output:
        model=train_results_dir / "model.json",
        predictions=train_results_dir / "predictions.csv",
    conda:
        envs_path("ebm.yml")
    log:
        train_log_dir / "decompose.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        python_path("decompose_model.py")


rule summarize_model:
    input:
        **rules.decompose_model.output,
        paths=rules.prepare_train_data.output.paths,
    output:
        train_results_dir / "summary.html",
    conda:
        envs_path("rmarkdown.yml")
    benchmark:
        train_results_dir / "summary.bench"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        rmd_path("train_summary.Rmd")


################################################################################
# prepare test data

prepare_x_file = "test_x.tsv"
prepare_log_file = "prepare.log"
prepare_bench_file = "prepare.bench"

test_dir = "test"
test_subdir = "{test_key}@{input_key}"

test_results_dir = train_results_dir / test_dir
test_log_dir = train_log_dir / test_dir

labeled_test_results_dir = test_results_dir / labeled_dir / test_subdir
unlabeled_test_results_dir = test_results_dir / unlabeled_dir / test_subdir

labeled_test_log_dir = test_log_dir / labeled_dir / test_subdir
unlabeled_test_log_dir = test_log_dir / unlabeled_dir / test_subdir


def test_data_input(annotated_path, wildcards):
    return {
        "annotated": expand(
            annotated_path,
            allow_missing=True,
            input_key=wildcards.test_key,
        ),
        "paths": rules.prepare_train_data.output.paths,
    }


# NOTE: the script will treat the input differently depending on if we want
# test_y (in this case it will perform all the label-related logic)
rule prepare_labeled_test_data:
    input:
        unpack(partial(test_data_input, rules.annotate_labeled_tsv.output)),
    output:
        test_x=labeled_test_results_dir / prepare_x_file,
        test_y=labeled_test_results_dir / "test_y.tsv",
    log:
        labeled_test_log_dir / prepare_log_file,
    benchmark:
        labeled_test_results_dir / prepare_bench_file
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        python_path("prepare_test.py")


use rule prepare_labeled_test_data as prepare_unlabeled_test_data with:
    input:
        unpack(partial(test_data_input, rules.annotate_unlabeled_tsv.output)),
    output:
        test_x=unlabeled_test_results_dir / prepare_x_file,
    log:
        unlabeled_test_log_dir / prepare_log_file,
    benchmark:
        unlabeled_test_results_dir / prepare_bench_file


################################################################################
# test model


test_bench_file = "test.bench"
test_log_file = "test.log"


def test_ebm_input(x_path):
    return {"model": rules.train_model.output.model, "test_x": x_path}


def test_ebm_output(test_path):
    return {
        "predictions": test_path / "predictions.csv",
        "explanations": test_path / "explanations.csv",
    }


rule test_labeled_ebm:
    input:
        **test_ebm_input(rules.prepare_labeled_test_data.output.test_x),
    output:
        **test_ebm_output(labeled_test_results_dir),
    conda:
        envs_path("ebm.yml")
    log:
        labeled_test_log_dir / test_log_file,
    resources:
        mem_mb=attempt_mem_gb(2),
    benchmark:
        labeled_test_results_dir / test_bench_file
    script:
        python_path("test_ebm.py")


use rule test_labeled_ebm as test_unlabeled_ebm with:
    input:
        **test_ebm_input(rules.prepare_unlabeled_test_data.output.test_x),
    output:
        **test_ebm_output(unlabeled_test_results_dir),
    log:
        unlabeled_test_log_dir / test_log_file,
    benchmark:
        unlabeled_test_results_dir / test_bench_file


################################################################################
# summarize test

test_summary_file = "summary.html"
test_summary_bench = "summary.bench"


rule summarize_labeled_test:
    input:
        **rules.test_labeled_ebm.output,
        truth_y=rules.prepare_labeled_test_data.output.test_y,
    output:
        labeled_test_results_dir / test_summary_file,
    conda:
        envs_path("rmarkdown.yml")
    benchmark:
        labeled_test_results_dir / test_summary_bench
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        rmd_path("test_summary.Rmd")


use rule summarize_labeled_test as summarize_unlabeled_test with:
    input:
        **rules.test_unlabeled_ebm.output,
    output:
        unlabeled_test_results_dir / test_summary_file,
    benchmark:
        unlabeled_test_results_dir / test_summary_bench


################################################################################
# global targets


def expand_unzip(path, keys, key_set):
    if isinstance(keys, dict):
        lookup = keys.values()
        keys = [*keys]
    else:
        lookup = keys
    subkeys = [[getattr(s, l) for l in lookup] for s in key_set]
    return expand(path, zip, **{k: [*s] for k, s in zip(keys, unzip(subkeys))})


def test_has_bench(run_keys_test):
    return (
        lookup_config(
            config,
            "inputs",
            run_keys_test.input_key,
            "test",
            run_keys_test.test_key,
            "benchmark",
        )
        is not None
    )


RunKeys = namedtuple("RunKeys", ["run_key", "filter_key", "input_keys", "inputs"])
RunKeysTrain = namedtuple(
    "RunKeysTrain",
    [
        "run_key",
        "filter_key",
        "input_keys",
        "input_key",
    ],
)
RunKeysTest = namedtuple(
    "RunKeysTest",
    [
        "run_key",
        "filter_key",
        "input_keys",
        "input_key",
        "test_key",
    ],
)


run_set = [
    RunKeys(k, f, INPUT_DELIM.join([*ns]), ns)
    for k, rs in config["ebm_runs"].items()
    for f in rs["filter"]
    for ns in rs["inputs"]
]

train_set = [
    RunKeysTrain(r.run_key, r.filter_key, r.input_keys, i)
    for r in run_set
    for i in r.inputs
]

test_set = [
    RunKeysTest(r.run_key, r.filter_key, r.input_keys, i, t)
    for r in run_set
    for i, ts in r.inputs.items()
    for t in ts
]
unlabeled_test_set, labeled_test_set = map(list, partition(test_has_bench, test_set))


def all_input_summary_files():
    def labeled_targets(key_set):
        return expand_unzip(
            rules.summarize_labeled_input.output,
            ["run_key", "filter_key", "input_key"],
            train_set,
        )

    train = labeled_targets(train_set)
    # labeled_test = labeled_targets(labeled_test_set)
    labeled_test = expand_unzip(
        rules.summarize_labeled_input.output,
        {
            "run_key": "run_key",
            "filter_key": "filter_key",
            "input_key": "test_key",
        },
        labeled_test_set,
    )
    unlabeled_test = expand_unzip(
        rules.summarize_unlabeled_input.output,
        {
            "run_key": "run_key",
            "filter_key": "filter_key",
            "input_key": "test_key",
        },
        unlabeled_test_set,
    )
    return train + labeled_test + unlabeled_test


def all_ebm_files():
    train = expand_unzip(
        rules.summarize_model.output,
        ["run_key", "filter_key", "input_keys"],
        run_set,
    )

    def expand_unzip_test(path, key_set):
        keys = ["run_key", "filter_key", "input_keys", "input_key", "test_key"]
        return expand_unzip(path, keys, key_set)

    # TODO these should eventually point to the test summary htmls
    labeled_test = expand_unzip_test(
        rules.summarize_labeled_test.output,
        labeled_test_set,
    )
    unlabeled_test = expand_unzip_test(
        rules.summarize_unlabeled_test.output,
        unlabeled_test_set,
    )

    return train + labeled_test + unlabeled_test


rule all_summary:
    input:
        all_input_summary_files(),


rule all_ebm:
    input:
        all_ebm_files(),
