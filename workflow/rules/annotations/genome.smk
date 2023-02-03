genome_dir = "genome"


# Download genome file (necessary from pretty much all the specific annotation
# rules). Use the first two columns of this table (chrom and length)
rule download_genome:
    output:
        annotations_src_dir / genome_dir / "chromInfo.txt.gz",
    params:
        url=partial(refkey_to_ref_wc, ["genome", "url"]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule get_genome:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_genome.output),
    output:
        annotations_tsv_dir / genome_dir / "genome.txt",
    log:
        annotations_log_dir / genome_dir / "genome.log",
    script:
        python_path("get_genome.py")
