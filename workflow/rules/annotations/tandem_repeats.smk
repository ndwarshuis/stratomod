tandem_repeats_src_dir = annotations_src_dir / "tandem_repeats"
tandem_repeats_results_dir = annotations_tsv_dir / "tandem_repeats"

# download this entire table as-is, we will select the right columns in a script
rule get_simreps_src:
    output:
        tandem_repeats_src_dir / "simple_repeats.tsv",
    params:
        url="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz",
    shell:
        "curl {params.url} | gunzip -c > {output}"

# NOTE sorting is done internally by the script
rule get_simple_reps:
    input:
        src=rules.get_simreps_src.output,
        genome=rules.get_genome.output,
    output:
        tandem_repeats_results_dir / "merged_simreps.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    script:
        str(scripts_dir / "get_simple_repeats.py") 
