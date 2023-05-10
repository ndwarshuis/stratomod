from scripts.python.common.config import attempt_mem_gb

map_dir = "mappability"
map_src = config.annotation_src_dir(log=False) / map_dir
map_src_log = config.annotation_src_dir(log=True) / map_dir
map_res = config.annotation_res_dir(log=False) / map_dir
map_res_log = config.annotation_res_dir(log=True) / map_dir


use rule download_labeled_query_vcf as download_mappability_high with:
    output:
        map_src / "high.bed.gz",
    log:
        map_src_log / "download_high.log",
    params:
        src=lambda w: config.references[w.ref_key].annotations.mappability.high.src,
    conda:
        "../../envs/utils.yml"
    localrule: True


use rule download_mappability_high as download_mappability_low with:
    output:
        map_src / "low.bed.gz",
    log:
        map_src_log / "download_low.log",
    params:
        src=lambda w: config.references[w.ref_key].annotations.mappability.low.src,
    localrule: True


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
        map_res_log / "mappability.log",
    log:
        map_res_log / "mappability.bench",
    conda:
        "../../envs/bio.yml"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        "../../scripts/python/bio/get_mappability_features.py"
