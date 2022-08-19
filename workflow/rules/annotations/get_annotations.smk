annotations_src_dir = resources_dir / "annotations"
annotations_tsv_dir = results_dir / "annotations"

# All these files from from here:
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

# Download genome file (necessary from pretty much all the specific annotation
# rules below). Use the first two columns of this table (chrom and length)
rule get_genome:
    output:
        annotations_src_dir / "genome.txt",
    params:
        url=config["resources"]["references"]["GRCh38"]["genome"]
    conda:
        envs_path("utils.yml")
    shell:
        f"""
        {sh_path('download_standardized')} gzip {{params.url}} 1 | \
        cut -f1,2 | \
        sort -k1,1n \
        > {{output}}
        """

include: "repeat_masker.smk"
include: "homopolymers.smk"
include: "mappability.smk"
include: "tandem_repeats.smk"
include: "segdups.smk"
