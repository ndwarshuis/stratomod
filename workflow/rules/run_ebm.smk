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
ebm_output_files = [
    ebm_dir / f
    for f in [
        "model.pickle",
        "train_x.pickle",
        "train_y.pickle",
        "test_x.pickle",
        "test_y.pickle",
        "config.yml",
    ]
]

################################################################################
# add annotations


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
        results_dir / "annotated_input" / "{input_key}" / "{filter_key}.tsv",
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
        ebm_output_files,
    params:
        config=lambda wildcards: json.dumps(lookup_ebm_run(wildcards)),
        out_dir=str(ebm_dir),
    conda:
        str(envs_dir / "ebm.yml")
    shell:
        """python workflow/scripts/run_ebm.py \
        -i {input} \
        -c '{params.config}' \
        -o {params.out_dir}
        """
