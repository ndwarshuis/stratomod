from scripts.python.common.config import attempt_mem_gb

mappability_dir = "mappability"
mappability_src_dir = annotations_src_dir / mappability_dir
mappability_results_dir = annotations_tsv_dir / mappability_dir
mappability_config_path = ["annotations", "mappability"]


rule download_mappability_high:
    output:
        mappability_src_dir / "high.bed.gz",
    params:
        url=lambda wildcards: config.refkey_to_annotations(
            wildcards.ref_key
        ).mappability.high.url,
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_mappability_high as download_mappability_low with:
    output:
        mappability_src_dir / "low.bed.gz",
    params:
        url=lambda wildcards: config.refkey_to_annotations(
            wildcards.ref_key
        ).mappability.low.url,


rule get_mappability:
    input:
        low=partial(
            expand_refkey_from_refsetkey,
            rules.download_mappability_low.output[0],
        ),
        high=partial(
            expand_refkey_from_refsetkey,
            rules.download_mappability_high.output[0],
        ),
    output:
        high=ensure(
            mappability_results_dir / "mappability_high.tsv.gz",
            non_empty=True,
        ),
        low=ensure(
            mappability_results_dir / "mappability_low_no_high.tsv.gz",
            non_empty=True,
        ),
    log:
        annotations_log_dir / mappability_dir / "mappability.log",
    conda:
        envs_path("bedtools.yml")
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        python_path("get_mappability_features.py")
