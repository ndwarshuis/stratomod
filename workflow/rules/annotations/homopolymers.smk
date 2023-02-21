from scripts.python.common.config import attempt_mem_gb

hp_dir = "homopolymers"
hp_res = config.annotation_dir(hp_dir, log=False)
hp_log = config.annotation_dir(hp_dir, log=True)


rule download_repseq:
    output:
        config.tool_resource_dir / "repseq.tar.gz",
    params:
        url=config.tools.repseq,
    conda:
        config.env_file("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(config.tool_dir(log=False) / "make" / "repseq"),
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        config.tool_dir(log=False) / "bin" / "repseq",
    conda:
        config.env_file("build")
    log:
        config.tool_dir(log=True) / "repseq.log",
    shell:
        "make -C {input} > {log} && mv {input}/repseq {output}"


# ASSUME the FASTA input to this is already standardized and filtered
rule find_simple_repeats:
    input:
        ref=partial(expand_refkey_from_refsetkey, rules.sdf_to_fasta.output),
        bin=rules.build_repseq.output,
    output:
        hp_res / "simple_repeats_p3.bed",
    benchmark:
        hp_log / "find_regions.bench"
    log:
        hp_log / "find_regions.log",
    resources:
        mem_mb=attempt_mem_gb(4),
    shell:
        """
        {input.bin} 1 4 {input.ref} 2> {log} | \
        sed '/^#/d' | \
        sort -k 1,1n -k 2,2n -k 3,3n \
        > {output}
        """


def homopolymer_file(ext):
    return wildcard_format_ext("homopolymers_{}", ["base"], ext)


rule get_homopolymers:
    input:
        bed=rules.find_simple_repeats.output,
        genome=rules.get_genome.output,
    output:
        ensure(hp_res / homopolymer_file("tsv.gz"), non_empty=True),
    conda:
        config.env_file("bedtools")
    log:
        hp_log / homopolymer_file("log"),
    benchmark:
        hp_log / homopolymer_file("bench")
    resources:
        mem_mb=attempt_mem_gb(16),
    script:
        config.python_script("bedtools/get_homopoly_features.py")
