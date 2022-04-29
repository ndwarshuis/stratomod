mappability_src_dir = annotations_src_dir / "mappability"
mappability_results_dir = annotations_tsv_dir / "mappability"


# ASSUME these are already sorted numerically and filtered for complete
# chromosomes


# As this file is being downloaded, strip the header, add a new column filled
# with 1's, and add our own header to identify this new column. When intersected
# with the feature tsv file, the 1's will be a binary representation of a
# variant being in a given mappability region (and non-overlaps will be empty,
# which presumably will be filled with 0's)
rule get_mappability_high_src:
    output:
        mappability_src_dir / "mappability_high.tsv",
    params:
        url=config["resources"]["annotations"]["mappability"]["high"],
        feature_name="MAP_difficult_250bp",
    shell:
        """
        echo 'chrom\tstart\tend\t{params.feature_name}' > {output}
        curl -sS -L {params.url} | \
        gunzip -c | \
        sed 's/$/\t1/' | \
        tail -n+2 \
        >> {output}
        """


use rule get_mappability_high_src as get_mappability_low_src with:
    output:
        mappability_src_dir / "mappability_low.tsv",
    params:
        url=config["resources"]["annotations"]["mappability"]["low"],
        feature_name="MAP_difficult_100bp",


rule subtract_high_from_low_mappability:
    input:
        low=rules.get_mappability_low_src.output,
        high=rules.get_mappability_high_src.output,
    output:
        mappability_results_dir / "mappability_low_no_high.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    script:
        str(scripts_dir / "get_mappability_features.py")
