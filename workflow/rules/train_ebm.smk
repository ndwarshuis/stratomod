import re
import json
import subprocess as sp
from more_itertools import flatten, partition
from scripts.common.config import attempt_mem_gb


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


git_tag = get_git_tag()

train_dir = (
    results_dir / "ebm" / ("%s-{input_keys}-{filter_key}-{run_key,[^/]+}" % git_tag)
)
train_log_dir = train_dir / "log"

input_delim = "&"

test_dir = train_dir / "test" / "{test_key}@{input_key}"
test_log_dir = test_dir / "log"

################################################################################
# add annotations


def annotation_input(tsv_path):
    return {
        "variants": tsv_path,
        "annotations": [
            rules.get_repeat_masker_classes.output,
            rules.get_tandem_repeats.output,
            rules.get_mappability_high_src.output,
            rules.subtract_high_from_low_mappability.output,
            rules.get_segdups.output,
            expand(
                rules.get_homopolymers.output,
                bases=config["features"]["homopolymers"]["bases"],
                allow_missing=True,
            ),
        ],
    }


annotated_dir = results_dir / "annotated_input"
labeled_annotated_dir = annotated_dir / "labeled" / "{input_key}"
unlabeled_annotated_dir = annotated_dir / "unlabeled" / "{input_key}"
annotated_tsv = "{filter_key}.tsv"
annotated_log = "{filter_key}.log"
annotated_bench = "{filter_key}.bench"


rule annotate_labeled_tsv:
    input:
        **annotation_input(rules.concat_tsv_files.output),
    output:
        labeled_annotated_dir / annotated_tsv,
    conda:
        str(envs_dir / "bedtools.yml")
    log:
        labeled_annotated_dir / annotated_log,
    benchmark:
        labeled_annotated_dir / annotated_bench
    resources:
        mem_mb=attempt_mem_gb(32),
    script:
        str(scripts_dir / "annotate.py")


use rule annotate_labeled_tsv as annotate_unlabeled_tsv with:
    input:
        **annotation_input(rules.parse_unlabeled_vcf.output),
    output:
        unlabeled_annotated_dir / annotated_tsv,
    log:
        unlabeled_annotated_dir / annotated_log,
    benchmark:
        unlabeled_annotated_dir / annotated_bench


################################################################################
# summarize annotated input

summary_output = "{filter_key}_summary.html"
summary_bench = "{filter_key}_summary.html"


rule summarize_labeled_input:
    input:
        rules.annotate_labeled_tsv.output,
    output:
        labeled_annotated_dir / summary_output,
    conda:
        str(envs_dir / "rmarkdown.yml")
    benchmark:
        labeled_annotated_dir / summary_bench
    resources:
        mem_mb=attempt_mem_gb(16),
    script:
        str(scripts_dir / "rmarkdown" / "input_summary.Rmd")


use rule summarize_labeled_input as summarize_unlabeled_input with:
    input:
        rules.annotate_unlabeled_tsv.output,
    output:
        unlabeled_annotated_dir / summary_output,
    benchmark:
        unlabeled_annotated_dir / summary_bench


def all_input_summary_files():
    input_keys, filter_keys = unzip(
        set(
            (i, f)
            for k, v in config["ebm_runs"].items()
            for f in v["filter"]
            for i in flatten(v["inputs"])
        )
    )
    return expand(
        [
            rules.summarize_labeled_input.output,
            rules.summarize_unlabeled_input.output,
        ],
        zip,
        input_key=[*input_keys],
        filter_key=[*filter_keys],
    )


rule all_summary:
    input:
        all_input_summary_files(),


################################################################################
# training
#
# assume that this will take care of test/train split, actual training, and
# pickling


rule prepare_train_data:
    input:
        lambda wildcards: expand(
            rules.annotate_labeled_tsv.output,
            allow_missing=True,
            input_key=wildcards.input_keys.split(input_delim),
        ),
    output:
        df=train_dir / "train.tsv",
        paths=train_dir / "input_paths.json",
    log:
        train_log_dir / "prepare.log",
    benchmark:
        train_log_dir / "prepare.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        str(scripts_dir / "prepare_train.py")


rule train_model:
    input:
        rules.prepare_train_data.output,
    output:
        **{
            n: str((train_dir / n).with_suffix(".csv"))
            for n in [
                "train_x",
                "train_y",
                "test_x",
                "test_y",
            ]
        },
        model=train_dir / "model.pickle",
        config=train_dir / "config.yml",
    conda:
        str(envs_dir / "ebm.yml")
    log:
        train_log_dir / "model.log",
    threads: 1
    # ASSUME total memory is proportional the number of EBMs that are trained in
    # parallel (see outer_bags parameter in the EBM function call)
    resources:
        mem_mb=lambda wildcards, threads, attempt: 16000 * threads * 2 ** (attempt - 1),
    benchmark:
        train_log_dir / "model.bench"
    script:
        str(scripts_dir / "train_ebm.py")


rule decompose_model:
    input:
        **rules.train_model.output,
    output:
        model=train_dir / "model.json",
        predictions=train_dir / "predictions.csv",
    conda:
        str(envs_dir / "ebm.yml")
    log:
        train_log_dir / "decompose.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        str(scripts_dir / "decompose_model.py")


rule summarize_model:
    input:
        **rules.decompose_model.output,
        paths=rules.prepare_train_data.output.paths,
    output:
        train_dir / "summary.html",
    conda:
        str(envs_dir / "rmarkdown.yml")
    benchmark:
        train_dir / "summary.bench"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        str(scripts_dir / "rmarkdown" / "model_summary.Rmd")


################################################################################
# prepare test data


def test_data_input(annotated_path, wildcards):
    return {
        "annotated": (
            lambda wildcards: expand(
                annotated_path,
                allow_missing=True,
                input_key=wildcards.test_key,
            ),
        ),
        "paths": rules.prepare_train_data.output.paths,
    }


# NOTE: the script will treat the input differently depending on if we want
# test_y (in this case it will perform all the label-related logic)
rule prepare_labeled_test_data:
    input:
        partial(test_data_input, rules.annotate_labeled_tsv.output),
    output:
        test_x=test_dir / "test_x.tsv",
        test_y=test_dir / "test_y.tsv",
    log:
        test_log_dir / "prepare.log",
    benchmark:
        test_log_dir / "prepare.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        str(scripts_dir / "prepare_test.py")


use rule prepare_labeled_test_data as prepare_unlabeled_test_data with:
    input:
        partial(test_data_input, rules.annotate_labeled_tsv.output),
    output:
        test_x=test_dir / "test_x.tsv",


################################################################################
# test model


def test_ebm_input(x_path):
    return {"model": rules.train_model.output.model, "test_x": test_x}


rule test_labeled_ebm:
    input:
        **test_ebm_input(rules.prepare_labeled_test_data.output.test_x),
    output:
        predictions=test_dir / "predictions.csv",
        explanations=test_dir / "explanation.csv",
    conda:
        str(envs_dir / "ebm.yml")
    log:
        test_log_dir / "model.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    benchmark:
        test_log_dir / "model.bench"
    script:
        str(scripts_dir / "test_ebm.py")


use rule test_labeled_ebm as test_unlabeled_ebm with:
    input:
        **test_ebm_input(rules.prepare_unlabeled_test_data.output.test_x),


# TODO add summary rmarkdown thingy here


################################################################################
# global targets


def all_ebm_files():
    def expand_unzip(path, wildcards, combinations):
        kwargs = {w: [*c] for w, c in zip(wildcards, unzip(combinations))}
        return expand(path, zip, **kwargs)

    def run_set_train(run_set):
        return expand_unzip(
            rules.summarize_model.output,
            ["run_key", "filter_key", "input_keys"],
            [(r, f, c) for r, f, _, c in run_set],
        )

    def test_has_bench(train_key, test_key):
        return config["inputs"][train_key]["test"][test_key]["benchmark"] is not None

    def expand_unzip_test(path, combinations):
        return expand_unzip(
            path,
            ["run_key", "filter_key", "input_keys", "input_key", "test_key"],
            combinations,
        )

    # TODO this should eventually point to the test summary htmls
    def run_set_test(run_set):
        test_set = [
            (k, f, c, i, t)
            for k, f, ns, c in run_set
            for n in ns
            for i, ts in ns.items()
            for t in ts
        ]
        labeled, unlabeled = partition(lambda s: test_has_bench(s[3], s[4]), test_set)
        return expand_unzip_test(
            rules.test_labeled_ebm.output, labeled
        ) + expand_unzip_test(rules.test_unlabeled_ebm.output, unlabeled)

    run_set = [
        (k, f, ns, input_delim.join([*ns]))
        for k, rs in config["ebm_runs"].items()
        for f in rs["filter"]
        for ns in rs["inputs"]
    ]
    return run_set_train(run_set) + run_set_test(run_set)


rule all_ebm:
    input:
        all_ebm_files(),
