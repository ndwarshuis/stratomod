from scripts.common.config import lookup_global_chr_filter

homopolymers_src_dir = annotations_src_dir / "homopolymers"
homopolymers_results_dir = annotations_tsv_dir / "homopolymers"
filtered_chrs = lookup_global_chr_filter(config)


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


# The main reason this rule is here is because I got tired of waiting for
# downstream steps to run for 30 minutes. If I filter to some small chromosome
# it makes testing waaaaay nicer.
rule filter_no_alt_analysis:
    input:
        rules.download_no_alt_analysis.output,
    output:
        homopolymers_results_dir / "GRCh38_no_alt_analysis_set_filtered.fa",
    conda:
        str(envs_dir / "biopython.yml")
    params:
        filt=filtered_chrs,
    script:
        str(scripts_dir / "filter_fasta.py")


def get_pasta():
    # all of it...
    n = (
        "filter_no_alt_analysis"
        if len(filtered_chrs) > 0
        else "download_no_alt_analysis"
    )
    return getattr(rules, n).output


rule find_simple_repeats:
    input:
        script=rules.download_find_regions_script.output,
        fasta=get_pasta(),
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


# This rule is here because I got tired of doing this step twice (once for AT
# and once for GC)
rule sort_and_filter_simple_repeats:
    input:
        rules.find_simple_repeats.output,
    output:
        homopolymers_results_dir / "simple_repeats_p3_sorted.bed",
    log:
        homopolymers_results_dir / "sorted.log",
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
    script:
        str(scripts_dir / "get_homopoly.py")
