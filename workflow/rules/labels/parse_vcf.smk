from functools import partial
from scripts.python.common.config import (
    lookup_config,
    attempt_mem_gb,
    flat_inputs,
    flat_chr_filters,
)

inputs_dir = resources_dir / "inputs"

input_results_dir = results_dir / "inputs" / "{input_key,[^/]+}"

prepare_dir = input_results_dir / "prepare"
labeled_dir = input_results_dir / "labeled"
unlabeled_dir = input_results_dir / "unlabeled"
alt_bench_dir = input_results_dir / "bench"

rtg_dir = labeled_dir / "rtg"

LABELS = ["fp", "fn", "tp"]


def lookup_input(wildcards, *args):
    return lookup_config(config, "inputs", wildcards.input_key, *args)


def lookup_train(wildcards, key):
    ns = flat_inputs(config)
    return ns[wildcards.input_key][key]


def lookup_chr_filter(wildcards):
    ref = lookup_train(wildcards, "ref")
    prefix = lookup_train(wildcards, "chr_prefix")
    fs = [
        ("X" if i == 23 else ("Y" if i == 24 else str(i)))
        for i in flat_chr_filters(config)[wildcards.input_key]
    ]
    # GRCh38 is like 'chr1, chr2, etc'
    pr = prefix if prefix is not None else ("chr" if ref == "GRCh38" else "")
    return "\\|".join([f"{pr}{f}\\b" for f in fs])


def lookup_train_bench(path, wildcards):
    return expand(path, bench_key=lookup_train(wildcards, "benchmark"))


include: "download_resources.smk"


################################################################################
# query vcf preprocessing


rule download_query_vcf:
    output:
        inputs_dir / "{input_key}.vcf",
    params:
        url=lambda wildcards: lookup_train(wildcards, "url"),
    conda:
        envs_path("utils.yml")
    shell:
        f"{sh_path('download_bedlike')} gzip {{params.url}} > {{output}}"


rule filter_query_vcf:
    input:
        rules.download_query_vcf.output,
    output:
        prepare_dir / "filtered.vcf",
    params:
        filt=lookup_chr_filter,
    conda:
        envs_path("utils.yml")
    shell:
        "sed -n '/^\(#\|{params.filt}\)/p' {input} > {output}"


# TODO this is (probably) just for DV VCFs
rule fix_refcall_query_vcf:
    input:
        rules.filter_query_vcf.output,
    output:
        prepare_dir / "fixed_refcall.vcf",
    conda:
        envs_path("utils.yml")
    shell:
        f"""
        cat {{input}} | \
        sed -e '/.RefCall./ s/\.\/\./0\/1/g' | \
        sed -e '/.RefCall./ s/0\/0/0\/1/g' \
        > {{output}}
        """


################################################################################
# vcf -> tsv (labeled)

# NOTE: Tabix won't work unless the input vcf is compressed, which is why we
# need to do this weird bgzip acrobatics with the query/bench


def rule_output_suffix(rule, suffix):
    return f"{getattr(rules, rule).output[0]}.{suffix}"


rule zip_query_vcf:
    input:
        rules.fix_refcall_query_vcf.output,
    output:
        rule_output_suffix("filter_query_vcf", "gz"),
    conda:
        envs_path("utils.yml")
    shell:
        "bgzip -c {input} > {output}"


rule generate_query_tbi:
    input:
        rules.zip_query_vcf.output,
    output:
        rule_output_suffix("zip_query_vcf", "tbi"),
    conda:
        envs_path("utils.yml")
    shell:
        "tabix -p vcf {input}"


use rule filter_query_vcf as filter_bench_vcf with:
    input:
        partial(lookup_train_bench, rules.download_bench_vcf.output),
    output:
        alt_bench_dir / "filtered.vcf",


use rule filter_query_vcf as filter_bench_bed with:
    input:
        partial(lookup_train_bench, rules.download_bench_bed.output),
    output:
        alt_bench_dir / "filtered.bed",


use rule zip_query_vcf as zip_bench_vcf with:
    input:
        rules.filter_bench_vcf.output,
    output:
        rule_output_suffix("filter_bench_vcf", "gz"),


use rule generate_query_tbi as generate_bench_tbi with:
    input:
        rules.zip_bench_vcf.output,
    output:
        rule_output_suffix("zip_bench_vcf", "tbi"),


# NOTE rtg won't write to a directory that already exists, so do this weird tmp
# file thing
# TODO add option to switch of the "--ref-overlap --all-records" thingy
# TODO --all-records = use all records including those that fail, make an option
# for this


rule label_vcf:
    input:
        query_vcf=rules.zip_query_vcf.output,
        query_tbi=rules.generate_query_tbi.output,
        bench_vcf=rules.zip_bench_vcf.output,
        bench_bed=rules.filter_bench_bed.output,
        bench_tbi=rules.generate_bench_tbi.output,
        sdf=lambda wildcards: expand(
            rules.download_ref_sdf.output,
            ref_key=lookup_train(wildcards, "ref"),
        ),
    output:
        [rtg_dir / f"{lbl}.vcf" for lbl in LABELS],
    conda:
        envs_path("rtg.yml")
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
        --no-gzip \
        -b {input.bench_vcf} \
        -e {input.bench_bed} \
        -c {input.query_vcf} \
        -o {params.tmp_dir} \
        -t {input.sdf} > {log} 2>&1

        mv {params.tmp_dir}/* {params.output_dir}

        rm -r {params.tmp_dir}
        """


rule parse_labeled_vcf:
    input:
        rtg_dir / "{label}.vcf",
    output:
        labeled_dir / "{filter_key}_{label}.tsv",
    log:
        labeled_dir / "{filter_key}_{label}.log",
    benchmark:
        labeled_dir / "{filter_key}_{label}.bench"
    conda:
        envs_path("bedtools.yml")
    script:
        python_path("parse_vcf_to_bed_ebm.py")


rule concat_labeled_tsvs:
    input:
        expand(rules.parse_labeled_vcf.output, label=LABELS, allow_missing=True),
    output:
        labeled_dir / "{filter_key}_labeled.tsv",
    conda:
        envs_path("bedtools.yml")
    benchmark:
        labeled_dir / "{filter_key}_concat.bench"
    resources:
        mem_mb=attempt_mem_gb(4),
    script:
        python_path("concat_tsv.py")


################################################################################
# vcf -> tsv (unlabeled)


use rule parse_labeled_vcf as parse_unlabeled_vcf with:
    input:
        rules.filter_query_vcf.output,
    output:
        unlabeled_dir / "{input_key}%{filter_key}.tsv",
    log:
        unlabeled_dir / "{input_key}%{filter_key}.log",
    benchmark:
        unlabeled_dir / "{input_key}%{filter_key}.bench"
