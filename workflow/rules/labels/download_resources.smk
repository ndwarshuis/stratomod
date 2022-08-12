from os.path import dirname
from functools import partial
from scripts.common.config import lookup_global_chr_filter

bench_dir = resources_dir / "bench"
ref_resources_dir = resources_dir / "reference"
ref_results_dir = results_dir / "reference"


def lookup_resource(*args):
    return lookup_config(config, "resources", *args)


def lookup_reference(wildcards):
    return lookup_resource("references", wildcards.ref_key, "sdf")


def lookup_benchmark(key, wildcards):
    return lookup_resource("benchmarks", wildcards.bench_key, key)


################################################################################
# download reference


rule download_ref_sdf:
    output:
        directory(ref_resources_dir / "{ref_key}.sdf"),
    params:
        url=lookup_reference,
        dir=lambda _, output: dirname(output[0]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -Ss {params.url} | bsdtar -xf - -C {params.dir}"


# TODO this isn't really downloading anything but couldn't think of where else
# to put it
rule sdf_to_fasta:
    input:
        rules.download_ref_sdf.output,
    output:
        ref_results_dir / "{ref_key}.fa",
    # if filter is empty, this will produce a blank string and rtg sdf2fasta
    # will filter nothing
    params:
        filt=" ".join(lookup_global_chr_filter(config)),
    conda:
        envs_path("rtg.yml")
    benchmark:
        ref_results_dir / "{ref_key}.bench",
    shell:
        """
        rtg sdf2fasta \
        -Z --line-length=70 -n \
        -i {input} \
        -o {output} {params.filt}
        """


################################################################################
# download benchmark files


def download_bench_vcf_cmd(wildcards, output):
    # dirty hack to fix the v4.2.1 benchmark (this should filter out the MHC
    # region on chr6, which we don't really want anyway)
    cmd = (
        "curl -Ss {u} | gunzip -c | grep -v 'GT:AD:PS' | bgzip -c > {o}"
        if wildcards.bench_key == "v4.2.1"
        else "curl -Ss -o {o} {u}"
    )
    return cmd.format(u=lookup_benchmark("vcf_url", wildcards), o=output)


rule download_bench_vcf:
    output:
        bench_dir / "{bench_key}.vcf.gz",
    params:
        cmd=download_bench_vcf_cmd,
    conda:
        envs_path("utils.yml")
    shell:
        "{params.cmd}"


rule download_bench_bed:
    output:
        bench_dir / "{bench_key}.bed",
    params:
        url=partial(lookup_benchmark, "bed_url"),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -Ss -o {output} {params.url}"


rule generate_bench_tbi:
    input:
        rules.download_bench_vcf.output,
    output:
        bench_dir / "{bench_key}.vcf.gz.tbi",
    conda:
        envs_path("utils.yml")
    shell:
        "tabix -p vcf {input}"
