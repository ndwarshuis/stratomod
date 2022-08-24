from scripts.python.common.config import (
    lookup_global_chr_filter,
    lookup_annotations,
    attempt_mem_gb,
)

tandem_repeats_dir = "tandem_repeats"
tandem_repeats_results_dir = annotations_tsv_dir / tandem_repeats_dir


# download this entire table as-is, we will select the right columns in a script
rule download_tandem_repeats:
    output:
        annotations_src_dir / tandem_repeats_dir / "simple_repeats.txt.gz",
    params:
        url=lookup_annotations(config)["simreps"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


# NOTE sorting is done internally by the script
rule get_tandem_repeats:
    input:
        src=rules.download_tandem_repeats.output,
        genome=rules.get_genome.output,
    output:
        tandem_repeats_results_dir / "tandem_repeats.tsv",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=lookup_global_chr_filter(config),
    log:
        annotations_log_dir / tandem_repeats_dir / "tandem_repeats.log",
    benchmark:
        tandem_repeats_results_dir / "tandem_repeats.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        python_path("get_tandem_repeat_features.py")
