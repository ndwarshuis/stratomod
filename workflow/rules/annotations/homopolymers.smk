from scripts.python.common.config import attempt_mem_gb

homopolymers_dir = "homopolymers"
homopolymers_src_dir = annotations_src_dir / homopolymers_dir
homopolymers_results_dir = annotations_tsv_dir / homopolymers_dir
homopolymers_log_dir = annotations_log_dir / homopolymers_dir


rule download_repseq:
    output:
        resources_dir / "tools" / "repseq.tar.gz",
    params:
        url=config["tools"]["repseq"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(results_dir / "tools" / "make" / "repseq"),
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        results_dir / "tools" / "bin" / "repseq",
    conda:
        envs_path("build.yml")
    log:
        log_dir / "tools" / "repseq.log",
    shell:
        "make -C {input} > {log} && mv {input}/repseq {output}"


# ASSUME the FASTA input to this is already standardized and filtered
rule find_simple_repeats:
    input:
        ref=partial(expand_refkey_from_refsetkey, rules.sdf_to_fasta.output),
        bin=rules.build_repseq.output,
    output:
        homopolymers_results_dir / "simple_repeats_p3.bed",
    conda:
        envs_path("find_simple_repeats.yml")
    benchmark:
        homopolymers_results_dir / "find_regions.bench"
    log:
        homopolymers_log_dir / "find_regions.log",
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
        homopolymers_results_dir / homopolymer_file("tsv.gz"),
    conda:
        envs_path("bedtools.yml")
    log:
        homopolymers_log_dir / homopolymer_file("log"),
    benchmark:
        homopolymers_results_dir / homopolymer_file("bench")
    resources:
        mem_mb=attempt_mem_gb(16),
    script:
        python_path("get_homopoly_features.py")
