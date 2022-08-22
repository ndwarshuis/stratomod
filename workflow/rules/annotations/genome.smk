genome_src_dir = annotations_src_dir / "genome"
genome_tsv_dir = annotations_tsv_dir / "genome"


# Download genome file (necessary from pretty much all the specific annotation
# rules). Use the first two columns of this table (chrom and length)
rule download_genome:
    output:
        genome_src_dir / "chromInfo.txt.gz",
    params:
        url=config["resources"]["references"]["GRCh38"]["genome"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_genome:
    input:
        rules.download_genome.output,
    output:
        genome_tsv_dir / "genome.txt",
    log:
        genome_tsv_dir / "genome.log",
    script:
        python_path("get_genome.py")
