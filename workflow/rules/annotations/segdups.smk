from scripts.common.config import lookup_global_chr_filter, lookup_annotations

segdups_src_dir = annotations_src_dir / "segdups"
segdups_results_dir = annotations_tsv_dir / "segdups"


# download this entire table as-is, we will select the right columns in a script
rule download_superdups:
    output:
        segdups_src_dir / "superdups.txt",
    params:
        url=lookup_annotations(config)["superdups"],
    shell:
        "curl -Ss {params.url} | gunzip -c > {output}"


# NOTE sorting is done internally by the script
rule get_segdups:
    input:
        rules.download_superdups.output,
    output:
        segdups_results_dir / "segdups.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    params:
        filt=lookup_global_chr_filter(config),
    log:
        segdups_results_dir / "segdups.log",
    script:
        str(scripts_dir / "get_segdup_features.py")
