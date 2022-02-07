import re
import json
import subprocess as sp


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


def lookup_ebm_run(wildcards):
    return config["ebm_runs"][wildcards.run_key]


git_tag = get_git_tag()

ebm_dir = results_dir / "ebm" / ("%s-{input_key}-{filter_key}-{run_key}" % git_tag)

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
        annotated_input_dir / "{filter_key}_summary.pdf",
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
            for i in v["inputs"]
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
        rules.add_annotations.output,
    output:
        ebm_dir / "input.tsv",
    params:
        config=lambda wildcards: json.dumps(lookup_ebm_run(wildcards)["features"]),
    shell:
        """
        python workflow/scripts/postprocess.py \
        -c '{params.config}' \
        -i {input} \
        -o {output}
        """


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
            n: str((ebm_dir / n).with_suffix(".pickle"))
            for n in [
                "model",
                "train_x",
                "train_y",
                "test_x",
                "test_y",
            ]
        },
        config=ebm_dir / "config.yml",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
        out_dir=str(ebm_dir),
    conda:
        str(envs_dir / "ebm.yml")
    script:
        str(scripts_dir / "run_ebm.py")


# shell:
#     """python workflow/scripts/run_ebm.py \
#     -i {input} \
#     -c '{params.config}' \
#     -o {params.out_dir}
#     """


rule summarize_ebm:
    input:
        **rules.train_ebm.output,
    output:
        ebm_dir / "model_summary.pdf",
    conda:
        str(envs_dir / "ebm.yml")
    script:
        str(scripts_dir / "rmarkdown" / "model_summary.Rmd")


def all_ebm_files():
    run_keys, input_keys, filter_keys = unzip(
        [
            (k, i, f)
            for k, v in config["ebm_runs"].items()
            for f in v["filter"]
            for i in v["inputs"]
        ]
    )
    return expand(
        rules.summarize_ebm.output,
        zip,
        run_key=[*run_keys],
        input_key=[*input_keys],
        filter_key=[*filter_keys],
    )


rule all_ebm:
    input:
        all_ebm_files(),
