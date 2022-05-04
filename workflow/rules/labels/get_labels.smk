from functools import partial
from scripts.common.config import lookup_config

inputs_dir = resources_dir / "inputs"

label_dir = results_dir / "labels" / "{input_key}"
alt_bench_dir = results_dir / "bench" / "{input_key}"
rtg_dir = label_dir / "rtg"

labels = ["fp", "fn", "tp"]


def lookup_input(wildcards, *args):
    return lookup_config(config, "inputs", wildcards.input_key, *args)


def lookup_chr_filter(wildcards):
    # make something that looks like 'chr1\b\|chr2\b...'
    return "\\|".join([f + "\\b" for f in lookup_input(wildcards, "chr_filter")])


def lookup_input_bench(path, wildcards):
    return expand(path, bench_key=lookup_input(wildcards, "benchmark"))


include: "download_resources.smk"


################################################################################
# VCF preprocessing

rule get_input_vcf:
    output:
        inputs_dir / "{input_key}.vcf.gz"
    params:
        url=lambda wildcards: lookup_input(wildcards, "url"),
    shell:
        "curl -o {output} {params.url}"


# TODO this is (probably) just for DV VCFs
rule preprocess_vcf:
    input:
        rules.get_input_vcf.output
    output:
        label_dir / "query.vcf.gz",
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        """
        gunzip -c {input} | \
        sed -e '/.RefCall./ s/\.\/\./0\/1/g' | \
        sed -e '/.RefCall./ s/0\/0/0\/1/g' | \
        bgzip -c > {output}
        """


rule index_vcf:
    input:
        rules.preprocess_vcf.output,
    output:
        label_dir / "query.vcf.gz.tbi",
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        "tabix -p vcf {input}"


rule filter_query_vcf:
    input:
        rules.preprocess_vcf.output,
    output:
        vcf=label_dir / "query_filtered.vcf.gz",
        tbi=label_dir / "query_filtered.vcf.gz.tbi",
    params:
        filt=lookup_chr_filter,
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        """
        gunzip -c {input} | \
        sed -n '/^\(#\|{params.filt}\)/p' | \
        bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}
        """


################################################################################
# bench preprocessing


use rule filter_query_vcf as filter_bench_vcf with:
    input:
        partial(lookup_input_bench, rules.get_bench_vcf.output),
    output:
        vcf=alt_bench_dir / "filtered.vcf.gz",
        tbi=alt_bench_dir / "filtered.vcf.gz.tbi",


rule filter_bench_bed:
    input:
        partial(lookup_input_bench, rules.get_bench_bed.output),
    output:
        alt_bench_dir / "filtered.bed",
    params:
        filt=lookup_chr_filter,
    shell:
        """
        sed -n '/^\(#\|{params.filt}\)/p' {input} > {output}
        """


################################################################################
# VCF -> tsv


# rtg won't output to a directory that already exists, so do this weird temp
# file thing
# TODO add option to switch of the "--ref-overlap --all-records" thingy
# TODO --all-records = use all records including those that fail, make an option
# for this


def get_truth_inputs(wildcards):
    # NOTE: tbi file not used in vcfeval CLI but still needed
    paths = (
        [
            rules.get_bench_vcf.output,
            rules.get_bench_bed.output,
            rules.get_bench_tbi.output,
        ]
        if len(lookup_input(wildcards, "chr_filter")) == 0
        else [
            rules.filter_bench_vcf.output.vcf,
            rules.filter_bench_bed.output,
            rules.filter_bench_vcf.output.tbi,
        ]
    )
    return {
        key: expand(
            path,
            input_key=wildcards.input_key,
            bench_key=lookup_input(wildcards, "benchmark"),
        )
        for key, path in zip(["truth_vcf", "truth_bed", "truth_tbi"], paths)
    }


def get_query_inputs(wildcards):
    # NOTE: tbi file not used in vcfeval CLI but still needed
    paths = (
        [
            rules.preprocess_vcf.output,
            rules.index_vcf.output,
        ]
        if len(lookup_input(wildcards, "chr_filter")) == 0
        else [
            rules.filter_query_vcf.output.vcf,
            rules.filter_query_vcf.output.tbi,
        ]
    )
    return dict(zip(["query_vcf", "query_bed"], paths))


rule get_vcf_labels:
    input:
        unpack(get_truth_inputs),
        unpack(get_query_inputs),
        sdf=lambda wildcards: expand(
            rules.get_ref_sdf.output,
            ref_key=lookup_input(wildcards, "ref"),
        ),
    output:
        [rtg_dir / f"{lbl}.vcf.gz" for lbl in labels],
    conda:
        str(envs_dir / "rtg.yml")
    params:
        extra="--ref-overlap --all-records",
        tmp_dir="/tmp/vcfeval",
        output_dir=lambda _, output: Path(output[0]).parent,
    log:
        rtg_dir / "vcfeval.log",
    shell:
        """
        rm -rf {params.tmp_dir}

        rtg vcfeval {params.extra} \
        -b {input.truth_vcf} \
        -e {input.truth_bed} \
        -c {input.query_vcf} \
        -o {params.tmp_dir} \
        -t {input.sdf} > {log} 2>&1

        mv {params.tmp_dir}/* {params.output_dir}

        rm -r {params.tmp_dir}
        """


rule unzip_vcf_labels:
    input:
        rtg_dir / "{label}.vcf.gz",
    output:
        label_dir / "{label}.vcf",
    shell:
        "gunzip {input} -c > {output}"


rule parse_label_vcf:
    input:
        rules.unzip_vcf_labels.output,
    output:
        label_dir / "{filter_key}_{label}.tsv",
    log:
        label_dir / "{filter_key}_{label}.log",
    shell:
        """
        python \
        workflow/scripts/parse_vcf_to_bed_ebm.py \
        --type {wildcards.filter_key} \
        --label {wildcards.label} \
        --input {input} \
        --output {output} > {log} 2>&1
        """


rule concat_tsv_files:
    input:
        expand(rules.parse_label_vcf.output, label=labels, allow_missing=True),
    output:
        label_dir / "{filter_key}_labeled.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    script:
        str(scripts_dir / "concat_tsv.py")


## TODO add filtering rules here if we wish
