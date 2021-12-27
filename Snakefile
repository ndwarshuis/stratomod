import pandas as pd
import json
from functools import partial
import subprocess as sp
from os.path import join
from snakemake.utils import min_version, validate

# include: "rules/blablabla.smk"

min_version("6.12")


configfile: "config/config.yml"


validate(config, "config/config-schema.yml")


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sp.run(args, capture_output=True).stdout.strip().decode()


def merge_dicts(d1, d2):
    # ASSUME: d2 is a subset of d1
    # NOTE: this will only recursively merge nested dicts (not lists, which
    # could also be present in a yaml config)
    def use_other_maybe(k, v):
        if k in d2:
            _v = d2[k]
            if isinstance(_v, dict):
                return merge_dicts(d1[k], _v)
            else:
                return _v
        else:
            return v

    if d1 is None:
        return d2
    if d2 is None:
        return d1
    return {k: use_other_maybe(k, v) for k, v in d1.items()}


def lookup_run_config(config, run_key):
    r = config["runs"][run_key]
    g = config["global"]
    return merge_dicts(g, r)


def lookup_run_json(config, run_key):
    return json.dumps(lookup_run_config(config, run_key))


def lookup_config(config, *keys):
    k = keys[0]
    ks = keys[1:]
    return config[k] if len(ks) == 0 else lookup_config(config[k], *ks)


git_tag = get_git_tag()
run_keys = list(config["runs"])

################################################################################
# output paths

resources_dir = "resources"

bench_dir = join(resources_dir, "bench")

annotations_dir = join(resources_dir, "annotations")

genome_file = join(annotations_dir, "genome.txt")

simple_repeat_file = join(annotations_dir, "simrep.tsv")
repeat_masker_file = join(annotations_dir, "rmsk.tsv")
homopolymer_at_file = join(annotations_dir, "homopoly_AT.tsv")
homopolymer_gc_file = join(annotations_dir, "homopoly_GC.tsv")
segdup_file = join(annotations_dir, "segdups.tsv")

results_dir = "results"

label_dir = join(results_dir, "labels")
preprocessing_dir = join(results_dir, "preprocessing")

ebm_dir = join(resources_dir, "ebm", "{}_{{run_key}}".format(git_tag))
ebm_input_tsv = join(ebm_dir, "input.tsv")
ebm_output_files = [
    join(ebm_dir, f)
    for f in [
        "model.pickle",
        "train_x.pickle",
        "train_y.pickle",
        "test_x.pickle",
        "text_y.pickle",
        "config.yml",
    ]
]

################################################################################
# main target
#
# Define what EBM models we want to run, and pin the output dirs to 'all'


rule all:
    input:
        expand(ebm_output_files, tag=git_tag, run_key=run_keys),


################################################################################
# VCF preprocessing


rule get_vcf:
    output:
        join(resources_dir, "query.vcf"),
    params:
        url=lookup_config(config, "resources", "query_url"),
    shell:
        "curl -O {params.url} | gunzip -c > {output}"


rule preprocess_vcf:
    input:
        rules.get_vcf.output,
    output:
        join(label_dir, "query_corrected_refcall.vcf"),
    shell:
        """
        cat {input} | \
        sed -e '/.RefCall./ s/\.\/\./0\/1/g' | \
        sed -e '/.RefCall./ s/0\/0/0\/1/g' \
        > {output}
        """


################################################################################
# VCF -> tsv


# TODO use more unique names for these file so it is less likely in the future
# that we will change the config and have an old/wrong file that won't be
# overwritten by snakemake
rule get_ref_sdf:
    output:
        join(resources_dir, "ref.sdf"),
    params:
        url=lookup_config(config, "resources", "ref", "sdf_url"),
    shell:
        "curl -O {params.url} | unzip -c > {output}"


rule get_bench_vcf:
    output:
        join(resources_dir, "bench.vcf"),
    params:
        url=lookup_config(config, "resources", "bench", "vcf_url"),
    shell:
        "curl -O {params.url} | gunzip -c > {output}"


rule get_bench_bed:
    output:
        join(resources_dir, "bench.bed"),
    params:
        url=lookup_config(config, "resources", "bench", "bed_url"),
    shell:
        "curl -o {output} {params.url}"


rule get_vcf_labels:
    input:
        query_vcf=rules.preprocess_vcf.output,
        truth_vcf=rules.get_bench_vcf.output,
        truth_bed=rules.get_bench_bed.output,
        sdf=rules.get_ref_sdf.output,
    output:
        tp="tp.vcf",
        fp="fp.vcf",
    params:
        extra="--ref-overlap --all-records",
    shell:
        """
        rtg vcfeval {params.extra} \
            -b {input.truth_vcf} \
            -e {input.truth_bed} \
            -c {input.query_vcf} \
            -o {output} \
        -t {input.sdf}"""


rule parse_label_vcf:
    input:
        rules.get_vcf_labels.output,
    output:
        tp="tp.tsv",
        fp="fp.tsv",
    shell:
        """python parse_script {input} --output {output}"""


rule concat_tsv_files:
    input:
        tp=rules.get_vcf_labels.output.tp,
        fp=rules.get_vcf_labels.output.fp,
    output:
        "labelled.tsv",
    shell:
        "cat {input.tp} {input.fp} > {output}"


## TODO add filtering rules here if we wish

################################################################################
# add superdups annotations


rule add_superdups:
    input:
        variants=rules.concat_tsv_files.output,
    output:
        join(preprocessing_dir, "superdups.tsv"),
    shell:
        """cat {input.variants} > {output}"""


################################################################################
# add homopolymer annotations


rule add_homopolymers:
    input:
        variants=rules.add_superdups.output,
    output:
        join(preprocessing_dir, "added_homopolymers.tsv"),
    shell:
        """cat {input.variants} > {output}"""


################################################################################
# add simple repeat annotations


rule get_genome_tsv:
    output:
        genome_file,
    shell:
        """
        mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
        -A -P 3306 -D hg38 -N -B -e \
        'SELECT chrom,size FROM chromInfo 
        WHERE chrom REGEXP "chr[0-9]{{1,2}}$"
        ORDER BY CAST(REPLACE(chrom,'chr','') as INT);' \
        > {output}
        """


rule get_simreps_tsv:
    output:
        simple_repeat_file,
    shell:
        """
        mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
        -A -P 3306 -D hg38 -B -e \
        'SELECT * FROM simpleRepeat WHERE chrom REGEXP "chr[0-9]{{1,2}}";' \
        > {output}
        """


rule add_simple_reps:
    input:
        variants=rules.add_homopolymers.output,
        annotations=rules.get_simreps_tsv.output,
        # genome=get_genome_file,
    output:
        join(preprocessing_dir, "added_simreps.tsv"),
    shell:
        ## TODO need to update the script to actually use the non-hardcoded genome
        """
        cat {input.variants} | python add_simreps {input.annotations} > {output}
        """


################################################################################
# add repeat masker annotations


rule get_repeat_masker_tsv:
    output:
        repeat_masker_file,
    shell:
        """
        mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
        -A -P 3306 -D hg38 -B -e \
        'SELECT genoName,genoStart,genoEnd,repClass FROM rmsk
        WHERE genoName REGEXP "chr[0-9]{{1,2}}";' \
        > {output}
        """


rule add_repeat_masker:
    input:
        variants=rules.add_homopolymers.output,
        annotations=rules.get_repeat_masker_tsv.output,
    output:
        join(preprocessing_dir, "added_repeat_masker.tsv"),
    shell:
        """
        cat {input.variants} | python add_simreps {input.annotations} > {output}
        """


################################################################################
# add mapability annotations


# rule add_mappability:
#     input:
#         add_repeat_masker.output,
#     output:
#         join(preprocessing_dir, "added_mappability.tsv"),


################################################################################
# postprocess output


def get_postprocess_config(wildcards):
    c = lookup_run_config(config, wildcards.run_key)
    return json.dumps(c["features"])


rule postprocess_output:
    input:
        rules.add_repeat_masker.output,
    output:
        ebm_input_tsv,
    params:
        config=get_postprocess_config,
    shell:
        """python some_script -c {params.config} -i {input} > {output}"""


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
        config=lambda wildcards: lookup_run_json(config, wildcards.run_key),
        out_dir=ebm_dir,
    shell:
        "python run_ebm -c '{params.config}' -o {params.out_dir}"
