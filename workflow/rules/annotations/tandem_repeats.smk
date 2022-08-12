from scripts.common.config import (
    lookup_global_chr_filter,
    lookup_annotations,
    attempt_mem_gb,
)

tandem_repeats_src_dir = annotations_src_dir / "tandem_repeats"
tandem_repeats_results_dir = annotations_tsv_dir / "tandem_repeats"


# download this entire table as-is, we will select the right columns in a script
rule get_simreps_src:
    output:
        tandem_repeats_src_dir / "simple_repeats.tsv",
    params:
        url=lookup_annotations(config)["simreps"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -Ss {params.url} | gunzip -c > {output}"


# NOTE sorting is done internally by the script
rule get_tandem_repeats:
    input:
        src=rules.get_simreps_src.output,
        genome=rules.get_genome.output,
    output:
        tandem_repeats_results_dir / "tandem_repeats.tsv",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=lookup_global_chr_filter(config),
    log:
        tandem_repeats_results_dir / "tandem_repeats.log",
    benchmark:
        tandem_repeats_results_dir / "tandem_repeats.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        scripts_path("get_tandem_repeat_features.py")
