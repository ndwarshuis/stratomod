from scripts.python.common.config import attempt_mem_gb

segdups_dir = "segdups"
segdups_src = config.annotation_resource_dir(segdups_dir)
segdups_tsv = config.annotation_dir(segdups_dir, log=True)
segdups_log = config.annotation_dir(segdups_dir, log=False)


rule download_superdups:
    output:
        segdups_src / "superdups.txt.gz",
    params:
        url=lambda wildcards: config.refkey_to_annotations(
            wildcards.ref_key
        ).superdups.url,
    conda:
        config.env_file("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_segdups:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_superdups.output),
    output:
        ensure(segdups_tsv / "segdups.tsv.gz", non_empty=True),
    conda:
        config.env_file("bedtools")
    params:
        filt=refsetkey_to_chr_indices_wc,
    log:
        segdups_log / "segdups.log",
    benchmark:
        segdups_log / "segdups.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        config.python_script("bedtools/get_segdup_features.py")
