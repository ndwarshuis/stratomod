genome_dir = "genome"


## ASSUME the fasta is sorted by chromosome index
rule get_genome:
    input:
        # partial(expand_refkey_from_refsetkey, rules.download_genome.output),
        rules.sdf_to_fasta.output,
    output:
        ensure(
            config.annotation_dir(genome_dir, log=False) / "genome.txt",
            non_empty=True,
        ),
    conda:
        config.env_file("utils")
    log:
        config.annotation_dir(genome_dir, log=True) / "genome.log",
    shell:
        """
        samtools faidx {input} -o - 2> {log} | \
        cut -f1,2 > \
        {output}
        """
