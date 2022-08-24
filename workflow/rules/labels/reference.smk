from os.path import dirname
from scripts.python.common.config import lookup_global_chr_filter

ref_resources_dir = resources_dir / "reference"
ref_results_dir = results_dir / "reference"


def lookup_reference(wildcards):
    return lookup_config(
        config,
        "resources",
        "references",
        wildcards.ref_key,
        "sdf",
    )


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
