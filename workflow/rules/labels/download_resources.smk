from os.path import dirname
from functools import partial

bench_dir = resources_dir / "bench"
ref_dir = resources_dir / "reference"


def lookup_resource(*args):
    return lookup_config(config, "resources", *args)


def lookup_reference(wildcards):
    return lookup_resource("references", wildcards.ref_key, "sdf")


def lookup_benchmark(key, wildcards):
    return lookup_resource("benchmarks", wildcards.bench_key, key)


################################################################################
# get reference sdf


rule get_ref_sdf:
    output:
        directory(ref_dir / "{ref_key}.sdf"),
    params:
        url=lookup_reference,
        dir=lambda _, output: dirname(output[0]),
    shell:
        "curl {params.url} | bsdtar -xf - -C {params.dir}"


################################################################################
# get benchmark files


def download_bench_vcf_cmd(wildcards, output):
    # dirty hack to fix the v4.2.1 benchmark
    cmd = (
        "curl {u} | gunzip -c | sed -e 's/GT:AD:PS/GT:PS/' | bgzip -c > {o}"
        if wildcards.bench_key == "v4.2.1"
        else "curl -o {o} {u}"
    )
    return cmd.format(u=lookup_benchmark("vcf_url", wildcards), o=output)


rule get_bench_vcf:
    output:
        bench_dir / "{bench_key}.vcf.gz",
    params:
        cmd=download_bench_vcf_cmd,
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        "{params.cmd}"


rule get_bench_bed:
    output:
        bench_dir / "{bench_key}.bed",
    params:
        url=partial(lookup_benchmark, "bed_url"),
    shell:
        "curl -o {output} {params.url}"


rule get_bench_tbi:
    output:
        bench_dir / "{bench_key}.vcf.gz.tbi",
    params:
        url=partial(lookup_benchmark, "tbi_url"),
    shell:
        "curl -o {output} {params.url}"
