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

ebm_dir = results_dir / "ebm" / ("%s-{input_keys}-{filter_key}-{run_key}" % git_tag)
ebm_log_dir = ebm_dir / "log"

input_delim = "&"

################################################################################
# add annotations


annotated_input_dir = results_dir / "annotated_input" / "{input_key}"


rule add_annotations:
    input:
        variants=rules.concat_tsv_files.output,
        tsvs=[
            expand(
                rules.get_homopolymers.output,
                bases=config["features"]["homopolymers"]["bases"],
                allow_missing=True,
            ),
            rules.get_repeat_masker_classes.output,
            rules.get_tandem_repeats.output,
            rules.get_mappability_high_src.output,
            rules.subtract_high_from_low_mappability.output,
            rules.get_segdups.output,
        ],
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
# postprocess output


rule postprocess_output:
    input:
        lambda wildcards: expand(
            rules.add_annotations.output,
            allow_missing=True,
            input_key=wildcards.input_keys.split(input_delim),
        ),
    output:
        df=ebm_dir / "input.tsv",
        paths=ebm_dir / "input_paths.json",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
    log:
        ebm_log_dir / "postprocess.log",
    benchmark:
        ebm_log_dir / "postprocess.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    script:
        str(scripts_dir / "postprocess.py")


################################################################################
# run EBM
#
# assume that this will take care of test/train split, actual training, and
# pickling


rule train_ebm:
    input:
        rules.postprocess_output.output,
    output:
        **{
            n: str((ebm_dir / n).with_suffix(".csv"))
            for n in [
                "train_x",
                "train_y",
                "test_x",
                "test_y",
            ]
        },
        model=ebm_dir / "model.pickle",
        config=ebm_dir / "config.yml",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
    conda:
        str(envs_dir / "ebm.yml")
    log:
        ebm_log_dir / "model.log",
    threads: 1
    # ASSUME total memory is proportional the number of EBMs that are trained in
    # parallel (see outer_bags parameter in the EBM function call)
    resources:
        mem_mb=lambda wildcards, threads, attempt: 16000 * threads * 2 ** (attempt - 1),
    benchmark:
        ebm_log_dir / "model.bench"
    script:
        str(scripts_dir / "run_ebm.py")


rule decompose_ebm:
    input:
        **rules.train_ebm.output,
    output:
        model=ebm_dir / "model.json",
        predictions=ebm_dir / "predictions.csv",
    conda:
        str(envs_dir / "ebm.yml")
    log:
        ebm_log_dir / "decompose.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        str(scripts_dir / "decompose_model.py")


rule summarize_ebm:
    input:
        **rules.decompose_ebm.output,
        paths=rules.postprocess_output.output.paths,
    output:
        ebm_dir / "summary.html",
    conda:
        str(envs_dir / "rmarkdown.yml")
    benchmark:
        ebm_dir / "summary.bench"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        str(scripts_dir / "rmarkdown" / "model_summary.Rmd")


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
