from scripts.python.common.config import lookup_annotations, attempt_mem_gb

mappability_dir = "mappability"
mappability_src_dir = annotations_src_dir / mappability_dir
mappability_results_dir = annotations_tsv_dir / mappability_dir

mappability_config = lookup_annotations(config)["mappability"]


rule download_mappability_high:
    output:
        mappability_src_dir / "mappability_high.bed.gz",
    params:
        url=mappability_config["high"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_mappability_high as download_mappability_low with:
    output:
        mappability_src_dir / "mappability_low.bed.gz",
    params:
        url=mappability_config["low"],


rule get_mappability:
    input:
        low=rules.download_mappability_low.output[0],
        high=rules.download_mappability_high.output[0],
    output:
        high=mappability_results_dir / "mappability_high.tsv",
        low=mappability_results_dir / "mappability_low_no_high.tsv",
    log:
        annotations_log_dir / mappability_dir / "mappability.log",
    conda:
        envs_path("bedtools.yml")
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        python_path("get_mappability_features.py")
