map_dir = "mappability"
map_src = config.features_src_dir(log=False) / map_dir
map_src_log = config.features_src_dir(log=True) / map_dir
map_res = config.features_res_dir(log=False) / map_dir
map_res_log = config.features_res_dir(log=True) / map_dir


rule download_mappability_high:
    output:
        map_src / "high.bed.gz",
    log:
        map_src_log / "download_high.log",
    params:
        src=lambda w: config.references[w.ref_key].feature_data.mappability.high.src,
    conda:
        "../../envs/utils.yml"
    localrule: True
    script:
        "../../scripts/python/bio/download_bed_or_vcf.py"


use rule download_mappability_high as download_mappability_low with:
    output:
        map_src / "low.bed.gz",
    log:
        map_src_log / "download_low.log",
    params:
        src=lambda w: config.references[w.ref_key].feature_data.mappability.low.src,
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
    script:
        "../../scripts/python/bio/get_mappability_features.py"
