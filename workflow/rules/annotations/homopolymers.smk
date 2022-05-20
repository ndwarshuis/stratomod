from scripts.common.config import lookup_annotations

homopolymers_src_dir = annotations_src_dir / "homopolymers"
homopolymers_results_dir = annotations_tsv_dir / "homopolymers"


rule find_simple_repeats:
    # TODO don't hardcode GRCh38 (when applicable)
    input:
        expand(rules.sdf_to_fasta.output, ref_key="GRCh38"),
    output:
        homopolymers_results_dir / "simple_repeats_p3.bed",
    conda:
        str(envs_dir / "find_simple_repeats.yml")
    params:
        script=str(scripts_dir / "find_regions.py"),
    benchmark:
        homopolymers_results_dir / "find_regions.bench",
    shell:
        """
        python {params.script} \
        -p 3 -d 100000 -t 100000 -q 100000 \
        {input} {output}
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
        str(envs_dir / "bedtools.yml")
    shell:
        """
        cat {input} | \
        python workflow/scripts/sort_and_filter_bed.py -c "#" \
        2> {log} > {output}
        """


rule get_homopolymers:
    input:
        bed=rules.sort_and_filter_simple_repeats.output,
        genome=rules.get_genome.output,
    output:
        homopolymers_results_dir / "homopolymers_{bases}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    log:
        homopolymers_results_dir / "homopolymers_{bases}.log",
    benchmark:
        homopolymers_results_dir / "homopolymers_{bases}.bench",
    resources:
        mem_mb=32000
    script:
        str(scripts_dir / "get_homopoly_features.py")
