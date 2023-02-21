genome_dir = "genome"


# Download genome file (necessary from pretty much all the specific annotation
# rules). Use the first two columns of this table (chrom and length)
rule download_genome:
    output:
        config.annotation_resource_dir(genome_dir) / "chromInfo.txt.gz",
    params:
        url=lambda wildcards: config.references[wildcards.ref_key].genome.url,
    conda:
        config.env_file("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_genome:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_genome.output),
    output:
        ensure(
            config.annotation_dir(genome_dir, log=False) / "genome.txt",
            non_empty=True,
        ),
    conda:
        config.env_file("bedtools")
    log:
        config.annotation_dir(genome_dir, log=True) / "genome.log",
    script:
        config.python_script("bedtools/get_genome.py")
