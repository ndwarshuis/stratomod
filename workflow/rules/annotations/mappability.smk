from scripts.python.common.config import attempt_mem_gb

map_dir = "mappability"
map_src = config.annotation_resource_dir(map_dir)
map_res = config.annotation_dir(map_dir, log=False)
map_log = config.annotation_dir(map_dir, log=True)


def map_lookup(wildcards):
    return config.refkey_to_annotations(wildcards.ref_key).mappability


rule download_mappability_high:
    output:
        map_src / "high.bed.gz",
    params:
        url=lambda wildcards: map_lookup(wildcards).high.url,
    conda:
        config.env_file("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_mappability_high as download_mappability_low with:
    output:
        map_src / "low.bed.gz",
    params:
        url=lambda wildcards: map_lookup(wildcards).low.url,


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
        high=ensure(map_res / "mappability_high.tsv.gz", non_empty=True),
        low=ensure(map_res / "mappability_low_no_high.tsv.gz", non_empty=True),
    log:
        map_log / "mappability.log",
    log:
        map_log / "mappability.bench",
    conda:
        config.env_file("bedtools")
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        config.python_script("bedtools/get_mappability_features.py")
