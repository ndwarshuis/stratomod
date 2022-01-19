rmsk_src_dir = annotations_src_dir / "repeat_masker"
rmsk_results_dir = annotations_tsv_dir / "repeat_masker"

# download the genoName, genoStart, genoEnd, repClass columns for this table
rule get_repeat_masker_src:
    output:
        rmsk_src_dir / "repeat_masker.tsv",
    params:
        url="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
    shell:
        """
        curl {params.url} | \
        gunzip -c | \
        cut -f6,7,8,12 \
        > {output}
        """

rmsk_classes = ["SINE", "LINE", "LTR", "Satellite"]

# NOTE sorting/filtering chromosomes is done internally by this script
rule get_repeat_masker_classes:
    input:
        rules.get_repeat_masker_src.output,
    output:
        expand(
            rmsk_results_dir / "repeat_masker_{cls}.tsv",
            cls=rmsk_classes,
        ),
    conda:
        "../envs/bedtools.yml"
    params:
        outdir=lambda _, output: Path(output[0]).parent,
        classes=",".join(rmsk_classes),
    shell:
        """
        python workflow/scripts/get_rmsk_classes.py \
        -i {input} \
        -o {params.outdir} \
        -c {params.classes}
        """
