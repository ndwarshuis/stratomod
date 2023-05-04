from functools import partial
import scripts.python.common.config as cfg
from scripts.python.common.functional import flip, compose


def lookup_benchmark_vcf(wildcards):
    cor = (
        config.refsetkey_to_ref(wildcards.refset_key)
        .benchmarks[wildcards.bench_key]
        .corrections.strip_IPS
    )
    return rules.fix_HG005_bench_vcf.output if cor else rules.filter_bench_vcf.output


def expand_benchmark_path(path, wildcards):
    return expand(
        path,
        allow_missing=True,
        ref_key=config.refsetkey_to_refkey(wildcards.refset_key),
    )


def rule_output_suffix(rulename, suffix):
    return f"{getattr(rules, rulename).output[0]}.{suffix}"


################################################################################
# reference


rule download_ref_sdf:
    output:
        directory(config.ref_resource_dir / "sdf"),
    params:
        url=lambda wildcards: config.references[wildcards.ref_key].sdf.src.url,
    conda:
        config.env_file("utils")
    shell:
        """
        mkdir {output} && \
        curl -Ss {params.url} | \
        bsdtar -xf - \
        --directory {output} \
        --strip-components=1
        """


rule sdf_to_fasta:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_ref_sdf.output),
    output:
        config.refset_dir(log=False) / "standardized_ref.fa.gz",
    params:
        filt=lambda wildcards: " ".join(
            config.refsetkey_to_sdf_chr_filter(wildcards.refset_key)
        ),
        prefix=lambda wildcards: config.refsetkey_to_ref(
            wildcards["refset_key"]
        ).sdf.chr_prefix,
    conda:
        config.env_file("rtg")
    benchmark:
        config.refset_dir(log=True) / "ref_standardized.bench"
    shell:
        """
        rtg sdf2fasta -Z --line-length=70 -n -i {input} -o - {params.filt} | \
        sed 's/^>{params.prefix}/>/' | \
        sed 's/^>X/>23/' | \
        sed 's/^>Y/>24/' |
        bgzip -c > {output}
        """


rule fasta_to_sdf:
    input:
        rules.sdf_to_fasta.output,
    output:
        directory(config.refset_dir(log=False) / "standardized_sdf"),
    conda:
        config.env_file("rtg")
    benchmark:
        config.refset_dir(log=True) / "sdf_standardized.bench"
    log:
        config.refset_dir(log=True) / "sdf_standardized.log",
    shell:
        "rtg format -o {output} {input} 2>&1 > {log}"


################################################################################
# download stratifications

# so far the only use for the stratifications here is to remove MHC since some
# benchmarks don't have VAF/DP here and removing only from the benchmark would
# produce automatic FPs


# TODO don't download anything here...this strat is 1 line long
rule write_mhc_strat:
    output:
        config.ref_resource_dir / "strats" / "mhc.bed.gz",
    params:
        regions=lambda w: config.references[w.ref_key].strats.mhc,
    localrule: True
    conda:
        config.env_file("bedtools")
    script:
        config.python_script("bedtools/write_bed.py")


################################################################################
# query vcf


rule download_labeled_query_vcf:
    output:
        config.labeled_query_resource_dir / cfg.wildcard_ext("l_query_key", "vcf.gz"),
    params:
        src=lambda w: config.labeled_queries[w.l_query_key].src,
    conda:
        config.env_file("utils")
    localrule: True
    script:
        config.python_script("bedtools/get_file.py")


use rule download_labeled_query_vcf as download_unlabeled_query_vcf with:
    output:
        config.unlabeled_query_resource_dir / cfg.wildcard_ext("ul_query_key", "vcf.gz"),
    params:
        src=lambda w: config.unlabeled_queries[w.ul_query_key].src,


rule filter_labeled_query_vcf:
    input:
        rules.download_labeled_query_vcf.output,
    output:
        config.query_prepare_dir(log=False, labeled=True) / "filtered.vcf",
    params:
        chr_prefix=lambda wildcards: config.querykey_to_chr_prefix(
            wildcards.l_query_key
        ),
    conda:
        config.env_file("bedtools")
    script:
        config.python_script("bedtools/standardize_bed.py")


use rule filter_labeled_query_vcf as filter_unlabeled_query_vcf with:
    input:
        rules.download_unlabeled_query_vcf.output,
    output:
        config.query_prepare_dir(log=False, labeled=False) / "filtered.vcf",
    params:
        chr_prefix=lambda wildcards: config.querykey_to_chr_prefix(
            wildcards.ul_query_key
        ),


# TODO this is (probably) just for DV VCFs
rule fix_refcall_query_vcf:
    input:
        rules.filter_labeled_query_vcf.output,
    output:
        config.query_prepare_dir(log=False, labeled=True) / "fixed_refcall.vcf",
    conda:
        config.env_file("utils")
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


rule zip_labeled_query_vcf:
    input:
        rules.fix_refcall_query_vcf.output,
    output:
        rule_output_suffix("fix_refcall_query_vcf", "gz"),
    conda:
        config.env_file("utils")
    shell:
        "bgzip -c {input} > {output}"


rule generate_query_tbi:
    input:
        rules.zip_labeled_query_vcf.output,
    output:
        rule_output_suffix("zip_labeled_query_vcf", "tbi"),
    conda:
        config.env_file("utils")
    shell:
        "tabix -p vcf {input}"


################################################################################
# benchmark vcf


use rule download_labeled_query_vcf as download_bench_vcf with:
    output:
        config.bench_resource_dir / cfg.wildcard_ext("bench_key", "vcf.gz"),
    params:
        src=lambda w: config.references[w.ref_key].benchmarks[w.bench_key].vcf.src,
    localrule: True


use rule filter_labeled_query_vcf as filter_bench_vcf with:
    input:
        partial(expand_benchmark_path, rules.download_bench_vcf.output),
    output:
        config.bench_dir(log=False) / "filtered.vcf",
    params:
        chr_prefix=lambda wildcards: config.benchkey_to_bed_chr_prefix(
            wildcards.refset_key, wildcards.bench_key
        ),


# NOTE: this avoids an error caused by vcfeval where it will strip out any
# fields in the SAMPLE column that end in a dot, which in turn will result in a
# FORMAT/SAMPLE cardinality mismatch.
rule fix_HG005_bench_vcf:
    input:
        rules.filter_bench_vcf.output,
    output:
        config.bench_dir(log=False) / "HG005_fixed.vcf",
    conda:
        config.env_file("utils")
    shell:
        """
        cat {input} | \
        sed 's/:IPS\t/\t/' | \
        sed 's/:[^:]\+$//' \
        > {output}
        """


use rule zip_labeled_query_vcf as zip_bench_vcf with:
    input:
        lookup_benchmark_vcf,
    output:
        config.bench_dir(log=False) / "final_bench.vcf.gz",


use rule generate_query_tbi as generate_bench_tbi with:
    input:
        rules.zip_bench_vcf.output,
    output:
        rule_output_suffix("zip_bench_vcf", "tbi"),


################################################################################
# benchmark bed


use rule download_labeled_query_vcf as download_bench_bed with:
    output:
        config.bench_resource_dir / cfg.wildcard_ext("bench_key", "bed.gz"),
    params:
        src=lambda w: config.references[w.ref_key].benchmarks[w.bench_key].bed.src,
    localrule: True


rule filter_bench_bed:
    input:
        partial(expand_benchmark_path, rules.download_bench_bed.output),
    output:
        config.bench_dir(log=False) / "filtered.bed",
    params:
        chr_prefix=lambda wildcards: config.benchkey_to_bed_chr_prefix(
            wildcards.refset_key, wildcards.bench_key
        ),
    conda:
        config.env_file("bedtools")
    script:
        config.python_script("bedtools/standardize_bed.py")


rule subtract_mhc_bench_bed:
    input:
        bed=rules.filter_bench_bed.output,
        mhc=partial(expand_refkey_from_refsetkey, rules.write_mhc_strat.output),
    output:
        config.bench_dir(log=False) / "noMHC.bed",
    conda:
        config.env_file("bedtools")
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


def vcf_bench_targets(wildcards):
    return {
        k: expand(
            v,
            refset_key=config.querykey_to_refsetkey(wildcards.l_query_key),
            bench_key=config.querykey_to_benchkey(wildcards.l_query_key),
        )
        for k, v in [
            ("bench_vcf", rules.zip_bench_vcf.output),
            ("bench_bed", rules.subtract_mhc_bench_bed.output),
            ("bench_tbi", rules.generate_bench_tbi.output),
        ]
    }


rule label_vcf:
    input:
        unpack(vcf_bench_targets),
        query_vcf=rules.zip_labeled_query_vcf.output,
        query_tbi=rules.generate_query_tbi.output,
        sdf=lambda wildcards: expand(
            rules.fasta_to_sdf.output,
            refset_key=config.querykey_to_refsetkey(wildcards.l_query_key),
        ),
    output:
        [config.vcfeval_dir(log=False) / f"{lbl}.vcf.gz" for lbl in cfg.VCFLabel.all()],
    conda:
        config.env_file("rtg")
    params:
        extra="--ref-overlap --all-records",
        tmp_dir=lambda wildcards: f"/tmp/vcfeval_{wildcards.l_query_key}",
        output_dir=lambda _, output: Path(output[0]).parent,
    log:
        config.vcfeval_dir(log=True) / "vcfeval.log",
    benchmark:
        config.vcfeval_dir(log=True) / "vcfeval.bench"
    resources:
        mem_mb=1000,
    threads: 1
    shell:
        """
        rm -rf {params.tmp_dir} && \

        rtg RTG_MEM=$(({resources.mem_mb}*80/100))M \
        vcfeval {params.extra} \
        --threads={threads} \
        -b {input.bench_vcf} \
        -e {input.bench_bed} \
        -c {input.query_vcf} \
        -o {params.tmp_dir} \
        -t {input.sdf} > {log} 2>&1 && \

        mv {params.tmp_dir}/* {params.output_dir} && \

        rm -r {params.tmp_dir}
        """


def labeled_file(ext):
    return cfg.wildcard_format_ext(f"{{}}_{{}}", ["filter_key", "label"], ext)


rule parse_labeled_vcf:
    input:
        config.vcfeval_dir(log=False) / cfg.wildcard_ext("label", "vcf.gz"),
    output:
        config.query_parsed_dir(labeled=True, log=False) / labeled_file("tsv.gz"),
    log:
        config.query_parsed_dir(labeled=True, log=True) / labeled_file("log"),
    benchmark:
        config.query_parsed_dir(labeled=True, log=True) / labeled_file("bench")
    conda:
        config.env_file("bedtools")
    params:
        query_key=lambda wildcards: wildcards.l_query_key,
    resources:
        mem_mb=cfg.attempt_mem_gb(2),
    script:
        config.python_script("bedtools/parse_vcf_to_bed_ebm.py")


rule concat_labeled_tsvs:
    input:
        expand(
            rules.parse_labeled_vcf.output,
            label=cfg.VCFLabel.all(),
            allow_missing=True,
        ),
    output:
        ensure(
            config.query_parsed_dir(labeled=True, log=False)
            / cfg.wildcard_format("{}_labeled.tsv.gz", "filter_key"),
            non_empty=True,
        ),
    conda:
        config.env_file("bedtools")
    benchmark:
        config.query_parsed_dir(labeled=True, log=True) / cfg.wildcard_format(
            "{}_concat.bench", "filter_key"
        )
    resources:
        mem_mb=cfg.attempt_mem_gb(4),
    script:
        config.python_script("bedtools/concat_tsv.py")


################################################################################
# vcf -> tsv (unlabeled)


def unlabeled_file(ext):
    return cfg.wildcard_ext("filter_key", ext)


use rule parse_labeled_vcf as parse_unlabeled_vcf with:
    input:
        rules.filter_unlabeled_query_vcf.output,
    output:
        ensure(
            config.query_parsed_dir(labeled=False, log=False)
            / unlabeled_file("tsv.gz"),
            non_empty=True,
        ),
    log:
        config.query_parsed_dir(labeled=False, log=True) / unlabeled_file("log"),
    params:
        query_key=lambda wildcards: wildcards.ul_query_key,
    resources:
        mem_mb=cfg.attempt_mem_gb(2),
    benchmark:
        config.query_parsed_dir(labeled=False, log=True) / unlabeled_file("bench")
