from scripts.python.common.config import attempt_mem_gb

segdups_dir = "segdups"
segdups_src = config.annotation_resource_dir(segdups_dir)
segdups_tsv = config.annotation_dir(segdups_dir, log=True)
segdups_log = config.annotation_dir(segdups_dir, log=False)


use rule download_labeled_query_vcf as download_superdups with:
    output:
        segdups_src / "superdups.txt.gz",
    params:
        src=lambda w: config.refkey_to_annotations(w.ref_key).superdups.src,
    localrule: True


rule get_segdups:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_superdups.output),
    output:
        ensure(segdups_tsv / "segdups.tsv.gz", non_empty=True),
    conda:
        config.env_file("bedtools")
    log:
        segdups_log / "segdups.log",
    benchmark:
        segdups_log / "segdups.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        config.python_script("bedtools/get_segdup_features.py")
