import pandas as pd
import subprocess as sp
from os.path import join
from snakemake.utils import min_version, validate

# include: "rules/blablabla.smk"

min_version("6.12")


configfile: "config/config.yml"


# validate(config, "config/schema.yml")


def get_git_tag():
    args = ["git", "describe", "--tags", "--abbrev=0", "--always"]
    return sb.run(args, capture_output=True).stdout.strip()


git_tag = get_git_tag()
cfg_names = []

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

ebm_output_dir = join(resources_dir, "ebm", "{git_tag}_{cfg_name}")
ebm_output_files = [
    join(ebm_output_dir, f)
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
        expand(EBM_OUTPUT_FILES_full, tag=git_tag, cfg_name=cfg_names),


################################################################################
# VCF preprocessing


rule get_vcf:
    output:
        join(resources_dir, "query.vcf"),
    shell:
        "curl blablabla"


rule preprocess_vcf:
    input:
        get_vcf.output,
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
    shell:
        """curl blablabla"""


rule get_bench_vcf:
    output:
        join(resources_dir, "bench.vcf"),
    shell:
        """curl blablabla"""


rule get_bench_bed:
    output:
        join(resources_dir, "bench.bed"),
    shell:
        """curl blablabla"""


rule get_vcf_labels:
    input:
        query_vcf=preprocess_vcf.output,
        truth_vcf=get_bench_vcf.output,
        truth_bed=get_bench_bed.output,
        sdf=get_ref_sdf.output,
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
        get_vcf_labels.output,
    output:
        tp="tp.tsv",
        fp="fp.tsv",
    shell:
        """python parse_script {input} --output {output}"""


rule concat_tsv_files:
    input:
        tp=get_vcf_labels.output.tp,
        fp=get_vcf_labels.output.fp,
    output:
        "labelled.tsv",
    shell:
        "cat {input.tp} {input.fp} > {output}"


## TODO add filtering rules here if we wish

################################################################################
# add superdups annotations


rule add_superdups:
    input:
        variants=concat_tsv_files.output,
    output:
        join(preprocessing_dir, "superdups.tsv"),
    shell:
        """cat {input.variants} > {output}"""


################################################################################
# add homopolymer annotations


rule add_homopolymers:
    input:
        variants=add_superdups.output,
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
        WHERE chrom REGEXP "chr[0-9]{1,2}$"
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
        'SELECT * FROM simpleRepeat WHERE chrom REGEXP "chr[0-9]{1,2}";' \
        > {output}
        """


rule add_simple_reps:
    input:
        variants=add_homopolymers.output,
        annotations=get_simreps_tsv.output,
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
        WHERE genoName REGEXP "chr[0-9]{1,2}";' \
        > {output}
        """


rule add_repeat_masker:
    input:
        variants=add_homopolymers.output,
        annotations=get_repeat_masker_tsv.output,
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


rule postprocess_output:
    input:
        add_repeat_masker.output,
    output:
        join(preprocessing_dir, "postprocess.tsv"),
    shell:
        """python some_postprocess_script {input} >> {output}"""


################################################################################
# run EBM
#
# assume that this will take care of test/train split, actual training, and
# pickling


rule train_ebm:
    input:
        postprocess_output.output,
    output:
        ebm_output_files,
    shell:
        "python run_ebm"
