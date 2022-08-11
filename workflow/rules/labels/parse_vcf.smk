from functools import partial
from scripts.common.config import (
    lookup_config,
    attempt_mem_gb,
    flat_inputs,
    flat_chr_filters,
)

inputs_dir = resources_dir / "inputs"

labeled_dir = results_dir / "labeled" / "{input_key}"
rtg_dir = labeled_dir / "rtg"

unlabeled_dir = results_dir / "unlabeled"

alt_bench_dir = results_dir / "bench" / "{input_key}"

LABELS = ["fp", "fn", "tp"]


def lookup_input(wildcards, *args):
    return lookup_config(config, "inputs", wildcards.input_key, *args)


def lookup_train(wildcards, key):
    ns = flat_inputs(config)
    return ns[wildcards.input_key][key]


# return lookup_input(wildcards, "train", key)


def lookup_chr_filter(wildcards):
    fs = flat_chr_filters(config)[wildcards.input_key]
    # make something that looks like 'chr1\b\|chr2\b...'
    # return "\\|".join([f + "\\b" for f in lookup_input(wildcards, "chr_filter")])
    return "\\|".join([f + "\\b" for f in fs])


def lookup_train_bench(path, wildcards):
    return expand(path, bench_key=lookup_train(wildcards, "benchmark"))


include: "download_resources.smk"


################################################################################
# query vcf preprocessing


rule download_input_vcf:
    output:
        inputs_dir / "{input_key}.vcf.gz",
    params:
        url=lambda wildcards: lookup_train(wildcards, "url"),
    conda:
        str(envs_dir / "utils.yml")
    shell:
        "curl -sS -o {output} {params.url}"


# TODO this is (probably) just for DV VCFs
rule fix_refcall_gt_field:
    input:
        rules.download_input_vcf.output,
    output:
        labeled_dir / "query.vcf.gz",
    conda:
        str(envs_dir / "utils.yml")
    shell:
        """
        gunzip -c {input} | \
        sed -e '/.RefCall./ s/\.\/\./0\/1/g' | \
        sed -e '/.RefCall./ s/0\/0/0\/1/g' | \
        bgzip -c > {output}
        """


rule index_vcf:
    input:
        rules.fix_refcall_gt_field.output,
    output:
        labeled_dir / "query.vcf.gz.tbi",
    conda:
        str(envs_dir / "utils.yml")
    shell:
        "tabix -p vcf {input}"


rule filter_query_vcf:
    input:
        rules.preprocess_vcf.output,
    output:
        vcf=labeled_dir / "query_filtered.vcf.gz",
        tbi=labeled_dir / "query_filtered.vcf.gz.tbi",
    params:
        filt=lookup_chr_filter,
    conda:
        str(envs_dir / "utils.yml")
    shell:
        """
        gunzip -c {input} | \
        sed -n '/^\(#\|{params.filt}\)/p' | \
        bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}
        """


################################################################################
# bench vcf preprocessing


use rule filter_query_vcf as filter_bench_vcf with:
    input:
        partial(lookup_train_bench, rules.get_bench_vcf.output),
    output:
        vcf=alt_bench_dir / "filtered.vcf.gz",
        tbi=alt_bench_dir / "filtered.vcf.gz.tbi",


rule filter_bench_bed:
    input:
        partial(lookup_train_bench, rules.get_bench_bed.output),
    output:
        alt_bench_dir / "filtered.bed",
    params:
        filt=lookup_chr_filter,
    shell:
        """
        sed -n '/^\(#\|{params.filt}\)/p' {input} > {output}
        """


################################################################################
# vcf -> tsv (labeled)


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
        if lookup_chr_filter(wildcards) == ""
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
            bench_key=lookup_train(wildcards, "benchmark"),
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
        if lookup_chr_filter(wildcards) == ""
        else [
            rules.filter_query_vcf.output.vcf,
            rules.filter_query_vcf.output.tbi,
        ]
    )
    return dict(zip(["query_vcf", "query_tbi"], paths))


rule label_vcf:
    input:
        unpack(get_truth_inputs),
        unpack(get_query_inputs),
        sdf=lambda wildcards: expand(
            rules.get_ref_sdf.output,
            ref_key=lookup_train(wildcards, "ref"),
        ),
    output:
        [rtg_dir / f"{lbl}.vcf.gz" for lbl in LABELS],
    conda:
        str(envs_dir / "rtg.yml")
    params:
        extra="--ref-overlap --all-records",
        tmp_dir=lambda wildcards: f"/tmp/vcfeval_{wildcards.input_key}",
        output_dir=lambda _, output: Path(output[0]).parent,
    log:
        rtg_dir / "vcfeval.log",
    benchmark:
        rtg_dir / "vcfeval.bench"
    resources:
        mem_mb=1000,
    threads: 1
    shell:
        """
        rm -rf {params.tmp_dir}

        rtg RTG_MEM=$(({resources.mem_mb}*80/100))M \
        vcfeval {params.extra} \
        --threads={threads} \
        -b {input.truth_vcf} \
        -e {input.truth_bed} \
        -c {input.query_vcf} \
        -o {params.tmp_dir} \
        -t {input.sdf} > {log} 2>&1

        mv {params.tmp_dir}/* {params.output_dir}

        rm -r {params.tmp_dir}
        """


rule unzip_labeled_vcf:
    input:
        rtg_dir / "{label}.vcf.gz",
    output:
        labeled_dir / "{label}.vcf",
    shell:
        "gunzip {input} -c > {output}"


rule parse_labeled_vcf:
    input:
        rules.unzip_labeled_vcf.output,
    output:
        labeled_dir / "{filter_key}_{label}.tsv",
    log:
        labeled_dir / "{filter_key}_{label}.log",
    benchmark:
        labeled_dir / "{filter_key}_{label}.bench"
    script:
        str(scripts_dir / "parse_vcf_to_bed_ebm.py")


rule concat_labeled_tsvs:
    input:
        expand(rules.parse_labeled_vcf.output, label=LABELS, allow_missing=True),
    output:
        labeled_dir / "{filter_key}_labeled.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    benchmark:
        labeled_dir / "{filter_key}_concat.bench"
    resources:
        mem_mb=attempt_mem_gb(4),
    script:
        str(scripts_dir / "concat_tsv.py")


################################################################################
# vcf -> tsv (unlabeled)


def get_query_input(wildcards):
    return (
        rules.preprocess_vcf.output
        if lookup_chr_filter(wildcards) == ""
        else rules.filter_query_vcf.output.vcf
    )


rule unzip_vcf:
    input:
        get_query_input,
    output:
        unlabeled_dir / "{input_key}.vcf",
    shell:
        "gunzip {input} -c > {output}"


# TODO need to make this work without a label argument
rule parse_unlabeled_vcf:
    input:
        rules.unzip_vcf.output,
    output:
        unlabeled_dir / "{input_key}.tsv",
    log:
        unlabeled_dir / "{input_key}.log",
    benchmark:
        unlabeled_dir / "{input_key}.bench"
    script:
        str(scripts_dir / "parse_vcf_to_bed_ebm.py")


# TODO add sort step here
