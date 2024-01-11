from scripts.python.common.config import wildcard_format_ext

hp_dir = "homopolymers"
hp_res = config.features_res_dir(log=False) / hp_dir
hp_log = config.features_res_dir(log=True) / hp_dir


rule download_repseq:
    output:
        config.tool_src_dir(log=False) / "repseq.tar.gz",
    params:
        url=config.tools.repseq,
    conda:
        "../../envs/utils.yml"
    localrule: True
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(config.tool_res_dir(log=False) / "make" / "repseq"),
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        config.tool_res_dir(log=False) / "bin" / "repseq",
    conda:
        "../../envs/build.yml"
    log:
        config.tool_res_dir(log=True) / "repseq.log",
    shell:
        "make -C {input} > {log} && mv {input}/repseq {output}"


# ASSUME the FASTA input to this is already standardized/filtered/sorted
rule find_simple_repeats:
    input:
        ref=partial(expand_refkey_from_refsetkey, rules.filter_sort_ref.output),
        bin=rules.build_repseq.output,
    output:
        hp_res / "simple_repeats_p3.bed.gz",
    benchmark:
        hp_log / "find_regions.bench"
    log:
        hp_log / "find_regions.log",
    shell:
        """
        gunzip -c {input.ref} | \
        {input.bin} 1 4 - 2> {log} | \
        sed '/^#/d' | \
        gzip -c > {output}
        """


def homopolymer_file(ext):
    return wildcard_format_ext("homopolymers_{}", ["base"], ext)


rule get_homopolymers:
    input:
        bed=rules.find_simple_repeats.output,
        genome=rules.fasta_to_genome.output,
    output:
        ensure(hp_res / homopolymer_file("tsv.gz"), non_empty=True),
    conda:
        "../../envs/bio.yml"
    log:
        hp_log / homopolymer_file("log"),
    benchmark:
        hp_log / homopolymer_file("bench")
    script:
        "../../scripts/python/bio/get_homopoly_features.py"
