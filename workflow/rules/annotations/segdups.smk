from scripts.python.common.config import attempt_mem_gb

segdups_dir = "segdups"
segdups_log = config.annotation_res_dir(log=False) / segdups_dir


use rule download_mappability_high as download_superdups with:
    output:
        config.annotation_src_dir(log=False) / segdups_dir / "superdups.txt.gz",
    log:
        config.annotation_src_dir(log=True) / segdups_dir / "download.log",
    params:
        src=lambda w: config.refkey_to_annotations(w.ref_key).superdups.src,
    localrule: True


rule get_segdups:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_superdups.output),
    output:
        ensure(
            config.annotation_res_dir(log=False) / segdups_dir / "segdups.tsv.gz",
            non_empty=True,
        ),
    conda:
        "../../envs/bio.yml"
    log:
        segdups_log / "segdups.log",
    benchmark:
        segdups_log / "segdups.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        "../../scripts/python/bio/get_segdup_features.py"
