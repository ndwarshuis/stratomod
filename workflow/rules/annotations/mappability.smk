from scripts.common.config import (
    lookup_annotations,
    fmt_mappability_feature,
    attempt_mem_gb,
)

mappability_src_dir = annotations_src_dir / "mappability"
mappability_results_dir = annotations_tsv_dir / "mappability"

mappability_config = lookup_annotations(config)["mappability"]


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
        url=mappability_config["high"],
        feature_name=fmt_mappability_feature(config, "high"),
    conda:
        str(envs_dir / "util.yml")
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
        url=mappability_config["low"],
        feature_name=fmt_mappability_feature(config, "low"),


rule subtract_high_from_low_mappability:
    input:
        low=rules.get_mappability_low_src.output,
        high=rules.get_mappability_high_src.output,
    output:
        mappability_results_dir / "mappability_low_no_high.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        str(scripts_dir / "get_mappability_features.py")
