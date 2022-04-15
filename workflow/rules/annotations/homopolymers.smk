homopolymers_src_dir = annotations_src_dir / "homopolymers"
homopolymers_results_dir = annotations_tsv_dir / "homopolymers"


rule download_no_alt_analysis:
    output:
        homopolymers_src_dir / "GRCh38_no_alt_analysis_set.fa",
    params:
        url=config["resources"]["annotations"]["homopolymers"]["fasta"],
    shell:
        "curl {params.url} | gunzip -c > {output}"


rule download_find_regions_script:
    output:
        homopolymers_src_dir / "find_regions.py",
    params:
        url=config["resources"]["annotations"]["homopolymers"]["find_regions"],
    shell:
        "curl -o {output} {params.url}"


rule find_simple_repeats:
    input:
        script=rules.download_find_regions_script.output,
        fasta=rules.download_no_alt_analysis.output,
    output:
        homopolymers_results_dir / "simple_repeats_p3.bed",
    conda:
        str(envs_dir / "find_simple_repeats.yml")
    shell:
        """
        python {input.script} \
        -p 3 -d 100000 -t 100000 -q 100000 \
        {input.fasta} {output}
        """


rule get_homopolymers:
    input:
        bed=rules.find_simple_repeats.output,
        genome=rules.get_genome.output,
    output:
        homopolymers_results_dir / "homopolymers_{bases}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    script:
        str(scripts_dir / "get_homopoly.py")
