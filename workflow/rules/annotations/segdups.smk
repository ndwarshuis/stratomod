from scripts.python.common.config import attempt_mem_gb

segdups_dir = "segdups"
segdups_results_dir = annotations_tsv_dir / segdups_dir


rule download_superdups:
    output:
        annotations_src_dir / segdups_dir / "superdups.txt.gz",
    params:
        url=partial(refkey_to_ref_wc, ["annotations", "superdups"]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_segdups:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_superdups.output),
    output:
        segdups_results_dir / "segdups.tsv.gz",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=refsetkey_to_chr_indices_wc,
    log:
        annotations_log_dir / segdups_dir / "segdups.log",
    benchmark:
        segdups_results_dir / "segdups.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        python_path("get_segdup_features.py")
