from scripts.python.common.config import lookup_annotations, attempt_mem_gb

homopolymers_src_dir = annotations_src_dir / "homopolymers"
homopolymers_results_dir = annotations_tsv_dir / "homopolymers"


rule find_simple_repeats:
    # TODO don't hardcode GRCh38 (when applicable)
    input:
        expand(rules.sdf_to_fasta.output, ref_key="GRCh38"),
    output:
        homopolymers_results_dir / "simple_repeats_p3.bed",
    conda:
        envs_path("find_simple_repeats.yml")
    benchmark:
        homopolymers_results_dir / "find_regions.bench"
    resources:
        mem_mb=attempt_mem_gb(4),
    shell:
        f"""
        python {python_path("find_regions.py")} \
        -p 3 -d 100000 -t 100000 -q 100000 \
        {{input}} {{output}}
        """


# This rule is here because I got tired of doing this step twice (once for AT
# and once for GC)
rule sort_and_filter_simple_repeats:
    input:
        rules.find_simple_repeats.output,
    output:
        homopolymers_results_dir / "simple_repeats_p3_sorted.bed",
    log:
        homopolymers_results_dir / "sorted.log",
    conda:
        envs_path("bedtools.yml")
    benchmark:
        homopolymers_results_dir / "sorted.bench"
    resources:
        mem_mb=attempt_mem_gb(16),
    shell:
        f"""
        cat {{input}} | \
        {sh_path('standardize_chrs')} 1 | \
        python {python_path('sort_and_filter_bed.py')} -c "#" \
        2> {{log}} > {{output}}
        """


rule get_homopolymers:
    input:
        bed=rules.sort_and_filter_simple_repeats.output,
        genome=rules.get_genome.output,
    output:
        homopolymers_results_dir / "homopolymers_{bases}.tsv",
    conda:
        envs_path("bedtools.yml")
    log:
        homopolymers_results_dir / "homopolymers_{bases}.log",
    benchmark:
        homopolymers_results_dir / "homopolymers_{bases}.bench"
    resources:
        mem_mb=attempt_mem_gb(16),
    script:
        python_path("get_homopoly_features.py")
