from scripts.python.common.config import (
    lookup_global_chr_filter,
    lookup_annotations,
    attempt_mem_gb,
)

segdups_src_dir = annotations_src_dir / "segdups"
segdups_results_dir = annotations_tsv_dir / "segdups"


rule download_superdups:
    output:
        segdups_src_dir / "superdups.txt.gz",
    params:
        url=lookup_annotations(config)["superdups"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_segdups:
    input:
        rules.download_superdups.output,
    output:
        segdups_results_dir / "segdups.tsv",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=lookup_global_chr_filter(config),
    log:
        segdups_results_dir / "segdups.log",
    benchmark:
        segdups_results_dir / "segdups.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        python_path("get_segdup_features.py")
