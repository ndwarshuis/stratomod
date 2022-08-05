import re
import json
import subprocess as sp
from more_itertools import flatten
from scripts.common.config import attempt_mem_gb


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


def lookup_ebm_run(wildcards):
    return config["ebm_runs"][wildcards.run_key]


git_tag = get_git_tag()

train_dir = results_dir / "ebm" / ("%s-{input_keys}-{filter_key}-{run_key}" % git_tag)
train_log_dir = train_dir / "log"

input_delim = "&"

test_dir = train_dir / "test" / "{test_key}"
test_log_dir = test_dir / "log"

################################################################################
# add annotations


annotated_input_dir = results_dir / "annotated_input" / "{input_key}"


rule add_annotations:
    input:
        variants=rules.concat_tsv_files.output,
        annotations=annotation_tsvs,
    output:
        annotated_input_dir / "{filter_key}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    log:
        annotated_input_dir / "{filter_key}.log",
    benchmark:
        annotated_input_dir / "{filter_key}.bench"
    resources:
        mem_mb=attempt_mem_gb(32),
    script:
        str(scripts_dir / "annotate.py")


################################################################################
# make summary table


rule make_input_summary:
    input:
        rules.add_annotations.output,
    output:
        annotated_input_dir / "{filter_key}_summary.html",
    conda:
        str(envs_dir / "rmarkdown.yml")
    benchmark:
        annotated_input_dir / "{filter_key}_summary.bench"
    resources:
        mem_mb=attempt_mem_gb(16),
    script:
        str(scripts_dir / "rmarkdown" / "input_summary.Rmd")


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
        rules.make_input_summary.output,
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
            rules.add_annotations.output,
            allow_missing=True,
            input_key=wildcards.input_keys.split(input_delim),
        ),
    output:
        df=train_dir / "train.tsv",
        paths=train_dir / "input_paths.json",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
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
        rules.postprocess_output.output,
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
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
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
# testing


# TODO this somehow needs to be made aware of which input key the predict key
# is under so that the path column can be populated
rule prepare_test_data:
    input:
        annotated=lambda wildcards: expand(
            rules.add_annotations.output,
            allow_missing=True,
            input_key=wildcards.predict_key,
        ),
        paths=rules.postprocess_outpput.output.paths,
    output:
        test_x=test_dir / "test_x.tsv",
        test_y=test_dir / "test_y.tsv",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
    log:
        test_log_dir / "prepare.log",
    benchmark:
        test_log_dir / "prepare.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        str(scripts_dir / "prepare_test.py")


rule test_ebm:
    input:
        model=rules.train_ebm.output.model,
        test_x=rules.prepare_test_data.output.test_x,
    output:
        predictions=test_dir / "predictions.csv",
        explanations=test_dir / "explanation.csv",
    conda:
        str(envs_dir / "ebm.yml")
    # TODO will likely need this later to get the features for headers
    # params:
    #     config=lambda wildcards: lookup_ebm_run(wildcards),
    log:
        test_log_dir / "model.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    benchmark:
        test_log_dir / "model.bench"
    script:
        str(scripts_dir / "test_ebm.py")


################################################################################
# global targets


def all_ebm_files():
    run_keys, combined_input_keys, filter_keys = unzip(
        [
            (k, input_delim.join(i), f)
            for k, v in config["ebm_runs"].items()
            for f in v["filter"]
            for i in v["inputs"]
        ]
    )
    return expand(
        rules.summarize_ebm.output,
        zip,
        run_key=[*run_keys],
        input_keys=[*combined_input_keys],
        filter_key=[*filter_keys],
    )


rule all_ebm:
    input:
        all_ebm_files(),
