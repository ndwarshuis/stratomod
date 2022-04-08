import re
import json
import subprocess as sp
from more_itertools import flatten


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


def lookup_ebm_run(wildcards):
    return config["ebm_runs"][wildcards.run_key]


git_tag = get_git_tag()

ebm_dir = results_dir / "ebm" / ("%s-{input_keys}-{filter_key}-{run_key}" % git_tag)

input_delim = "&"

################################################################################
# add annotations


annotated_input_dir = results_dir / "annotated_input" / "{input_key}"


rule add_annotations:
    input:
        variants=rules.concat_tsv_files.output,
        tsvs=[
            rules.get_repeat_masker_classes.output,
            rules.get_simple_reps.output,
            rules.get_mappability_high_src.output,
            rules.get_mappability_low_src.output,
            expand(
                rules.get_segdups.output,
                colname=list(segdups_cols),
                allow_missing=True,
            ),
            expand(
                rules.get_homopolymers.output,
                bases=["AT", "GC"],
                allow_missing=True,
            ),
        ],
    output:
        annotated_input_dir / "{filter_key}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    shell:
        """
        python workflow/scripts/annotate.py \
        -i {input.variants} \
        -t {input.tsvs} \
        -o {output}
        """


################################################################################
# make summary table


rule make_input_summary:
    input:
        rules.add_annotations.output,
    output:
        annotated_input_dir / "{filter_key}_summary.html",
    conda:
        str(envs_dir / "rmarkdown.yml")
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
        paths=ebm_dir / "input_paths.yml",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
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
    script:
        str(scripts_dir / "decompose_model.py")


rule summarize_ebm:
    input:
        **rules.decompose_ebm.output,
    output:
        ebm_dir / "summary.html",
    conda:
        str(envs_dir / "rmarkdown.yml")
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
