from functools import partial
from scripts.python.common.config import (
    attempt_mem_gb,
    inputkey_to_bench_correction,
)
from scripts.python.common.functional import flip, compose


include: "reference.smk"


# resource dirs

inputs_dir = resources_dir / "inputs"
bench_dir = resources_dir / "bench" / all_wildcards["ref_key"]

# relative output dirs

rel_input_results_dir = (
    Path("inputs") / all_wildcards["refset_key"] / all_wildcards["input_key"]
)
rel_prepare_dir = rel_input_results_dir / "prepare"
rel_labeled_dir = rel_input_results_dir / "label"
rel_rtg_dir = rel_labeled_dir / "rtg"
rel_unlabeled_dir = rel_input_results_dir / "unlabeled"
rel_alt_bench_dir = rel_input_results_dir / "bench"

# results dirs

prepare_dir = results_dir / rel_prepare_dir
labeled_dir = results_dir / rel_labeled_dir
unlabeled_dir = results_dir / rel_unlabeled_dir
alt_bench_dir = results_dir / rel_alt_bench_dir
rtg_dir = results_dir / rel_rtg_dir

# log dirs


def lookup_benchmark_vcf(wildcards):
    return (
        rules.fix_HG005_bench_vcf.output
        if inputkey_to_bench_correction(config, "strip_IPS", wildcards.input_key)
        else rules.filter_bench_vcf.output
    )


# this is necessary because benchmark output files depend on both the input key
# and the refset key, so both need to be expanded when referring to the input
# path (which are in terms of the ref key and bench key)
def expand_benchmark_path(path, wildcards):
    return compose(
        partial(flip(expand_benchkey_from_inputkey), wildcards),
        partial(flip(expand_refkey_from_refsetkey), wildcards),
    )(path)


################################################################################
# query vcf


rule download_query_vcf:
    output:
        inputs_dir / wildcard_ext("input_key", "vcf.gz"),
    params:
        url=partial(inputkey_to_input_wc, ["url"]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule filter_query_vcf:
    input:
        rules.download_query_vcf.output,
    output:
        prepare_dir / "filtered.vcf",
    params:
        gzip_in=True,
        gzip_out=False,
    script:
        python_path("standardize_bed.py")


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
# query vcf label preprocessing

# NOTE: Tabix won't work unless the input vcf is compressed, which is why we
# need to do this weird bgzip acrobatics with the query/bench


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


################################################################################
# benchmark vcf


rule download_bench_vcf:
    output:
        bench_dir / wildcard_ext("bench_key", "vcf.gz"),
    params:
        url=partial(refkey_to_benchmark_wc, "vcf_url"),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule filter_query_vcf as filter_bench_vcf with:
    input:
        partial(expand_benchmark_path, rules.download_bench_vcf.output),
    output:
        alt_bench_dir / "filtered.vcf",


# NOTE: this avoids an error caused by vcfeval where it will strip out any
# fields in the SAMPLE column that end in a dot, which in turn will result in a
# FORMAT/SAMPLE cardinality mismatch.
rule fix_HG005_bench_vcf:
    input:
        rules.filter_bench_vcf.output,
    output:
        alt_bench_dir / "HG005_fixed.vcf",
    conda:
        envs_path("utils.yml")
    shell:
        """
        cat {input} | \
        sed 's/:IPS\t/\t/' | \
        sed 's/:[^:]\+$//' \
        > {output}
        """


use rule zip_query_vcf as zip_bench_vcf with:
    input:
        lookup_benchmark_vcf,
    output:
        alt_bench_dir / "final_bench.vcf.gz",


use rule generate_query_tbi as generate_bench_tbi with:
    input:
        rules.zip_bench_vcf.output,
    output:
        rule_output_suffix("zip_bench_vcf", "tbi"),


################################################################################
# benchmark bed


rule download_bench_bed:
    output:
        bench_dir / wildcard_ext("bench_key", "bed"),
    params:
        url=partial(refkey_to_benchmark_wc, "bed_url"),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule filter_bench_bed:
    input:
        partial(expand_benchmark_path, rules.download_bench_bed.output),
    output:
        alt_bench_dir / "filtered.bed",
    params:
        gzip_in=False,
        gzip_out=False,
    script:
        python_path("standardize_bed.py")


rule standardize_mhc_strat:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_mhc_strat.output),
    output:
        alt_bench_dir / "strats" / "mhc_standardized.bed.gz",
    params:
        gzip_in=True,
        gzip_out=True,
    script:
        python_path("standardize_bed.py")


rule subtract_mhc_bench_bed:
    input:
        bed=rules.filter_bench_bed.output,
        mhc=partial(expand_refkey_from_refsetkey, rules.standardize_mhc_strat.output),
    output:
        alt_bench_dir / "noMHC.bed",
    output:
        prepare_dir / "no_mhc.vcf",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip {input.mhc} -c | \
        bedtools subtract -a {input.bed} -b - \
        > {output}
        """


################################################################################
# vcf -> tsv (labeled)

# NOTE rtg won't write to a directory that already exists, so do this weird tmp
# file thing
# TODO add option to switch of the "--ref-overlap --all-records" thingy
# TODO --all-records = use all records including those that fail, make an option
# for this
# TODO rtg apparently outputs a log file on its own with everything else (which
# has more in it than piping stdout/stderr as below)


rule label_vcf:
    input:
        query_vcf=rules.zip_query_vcf.output,
        query_tbi=rules.generate_query_tbi.output,
        bench_vcf=rules.zip_bench_vcf.output,
        bench_bed=rules.subtract_mhc_bench_bed.output,
        bench_tbi=rules.generate_bench_tbi.output,
        sdf=partial(expand_refkey_from_refsetkey, rules.fasta_to_sdf.output),
    output:
        [rtg_dir / f"{lbl}.vcf.gz" for lbl in ALL_LABELS],
    conda:
        envs_path("rtg.yml")
    params:
        extra="--ref-overlap --all-records",
        tmp_dir=lambda wildcards: f"/tmp/vcfeval_{wildcards.input_key}",
        output_dir=lambda _, output: Path(output[0]).parent,
    log:
        log_dir / rel_rtg_dir / "vcfeval.log",
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
        -b {input.bench_vcf} \
        -e {input.bench_bed} \
        -c {input.query_vcf} \
        -o {params.tmp_dir} \
        -t {input.sdf} > {log} 2>&1

        mv {params.tmp_dir}/* {params.output_dir}

        rm -r {params.tmp_dir}
        """


def labeled_file(ext):
    return wildcard_format_ext(f"{{}}_{{}}", ["filter_key", "label"], ext)


rule parse_labeled_vcf:
    input:
        rtg_dir / wildcard_ext("label", "vcf.gz"),
    output:
        labeled_dir / labeled_file("tsv.gz"),
    log:
        log_dir / rel_labeled_dir / labeled_file("log"),
    benchmark:
        labeled_dir / labeled_file("bench")
    conda:
        envs_path("bedtools.yml")
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        python_path("parse_vcf_to_bed_ebm.py")


rule concat_labeled_tsvs:
    input:
        expand(rules.parse_labeled_vcf.output, label=ALL_LABELS, allow_missing=True),
    output:
        labeled_dir / wildcard_format("{}_labeled.tsv.gz", "filter_key"),
    conda:
        envs_path("bedtools.yml")
    benchmark:
        labeled_dir / wildcard_format("{}_concat.bench", "filter_key")
    resources:
        mem_mb=attempt_mem_gb(4),
    script:
        python_path("concat_tsv.py")


################################################################################
# vcf -> tsv (unlabeled)


def unlabeled_file(ext):
    return wildcard_format_ext("{}%{}", ["input_key", "filter_key"], ext)


use rule parse_labeled_vcf as parse_unlabeled_vcf with:
    input:
        rules.filter_query_vcf.output,
    output:
        unlabeled_dir / unlabeled_file("tsv.gz"),
    log:
        log_dir / rel_unlabeled_dir / unlabeled_file("log"),
    resources:
        mem_mb=attempt_mem_gb(2),
    benchmark:
        unlabeled_dir / unlabeled_file("bench")
