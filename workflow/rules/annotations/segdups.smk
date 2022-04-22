segdups_src_dir = annotations_src_dir / "segdups"
segdups_results_dir = annotations_tsv_dir / "segdups"

segdups_stats = ["min", "max", "count", "mean"]
segdups_cols = {"alignL": 19, "fracMatchIndel": 28}


rule download_superdups:
    output:
        segdups_src_dir / "superdups.txt",
    params:
        url=config["resources"]["annotations"]["superdups"],
    shell:
        "curl {params.url} | gunzip -c > {output}"


# TODO my current self does not understand why my past self wrote a shell
# command in this manner...my future self should fix that
rule get_segdups:
    input:
        rules.download_superdups.output,
    output:
        segdups_results_dir / "merged_segdups_{colname}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    params:
        stats=",".join(segdups_stats),
        col=lambda wildcards: segdups_cols[wildcards.colname],
        header=lambda wildcards: "\t".join(
            ["chr", "start", "end"]
            + [f"{wildcards.colname}_{s}" for s in segdups_stats]
        ),
    log:
        segdups_results_dir / "merged_segdups_{colname}.log",
    shell:
        """
        echo '{params.header}' > {output}

        cat {input} | \
        cut -f2-4,{params.col} | \
        python workflow/scripts/sort_and_filter_bed.py | \
        mergeBed -i stdin -c 4 -o {params.stats} \
        2> {log} >> {output}
        """
