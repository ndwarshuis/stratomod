from os.path import dirname
from functools import partial
from scripts.python.common.config import lookup_global_chr_filter

ref_resources_dir = resources_dir / "reference" / "{ref_key}"
ref_results_dir = results_dir / "reference" / "{ref_key}"


def lookup_reference_key(which, wildcards):
    return lookup_config(
        config,
        "resources",
        "references",
        wildcards.ref_key,
        which,
    )

def lookup_reference(wildcards):
    return lookup_reference_key("sdf", wildcards)


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


rule sdf_to_fasta:
    input:
        rules.download_ref_sdf.output,
    output:
        ref_results_dir / "{ref_key}.fa",
    # if filter is empty, this will produce a blank string and rtg sdf2fasta
    # will filter nothing
    params:
        # TODO this will need to change to accommodate other references
        filt=" ".join([f"chr{i}" for i in lookup_global_chr_filter(config)]),
    conda:
        envs_path("rtg.yml")
    benchmark:
        ref_results_dir / "{ref_key}.bench"
    shell:
        """
        rtg sdf2fasta \
        -Z --line-length=70 -n \
        -i {input} \
        -o {output} {params.filt}
        """

################################################################################
# download stratifications

# so far the only use for the stratifications here is to remove MHC since some
# benchmarks don't have VAF/DP here and removing only from the benchmark would
# produce automatic FPs


rule download_mhc_strat:
    output:
        ref_resources_dir / "strats" / "mhc.bed.gz"
    params:
        url=lambda wildcards: lookup_reference_key("strats", wildcards)["mhc"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"

