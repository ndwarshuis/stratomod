from scripts.common.config import lookup_global_chr_filter, lookup_annotations

rmsk_src_dir = annotations_src_dir / "repeat_masker"
rmsk_results_dir = annotations_tsv_dir / "repeat_masker"
rmsk_result_prefix = "repeat_masker"

# TODO move this to the config
rmsk_classes = {
    "SINE": [],
    "LINE": ["L1", "L2", "CR1", "RTE-X", "RTE-BovB", "Penelope", "Dong-R4"],
    "LTR": [],
    "Satellite": [],
}


# download the genoName, genoStart, genoEnd, repClass columns for this table
rule get_repeat_masker_src:
    output:
        rmsk_src_dir / "repeat_masker.tsv",
    params:
        url=lookup_annotations(config)["repeat_masker"]
    shell:
        "curl -Ss {params.url} | gunzip -c > {output}"


# NOTE sorting/filtering chromosomes is done internally by this script
rule get_repeat_masker_classes:
    input:
        rules.get_repeat_masker_src.output,
    output:
        [rmsk_results_dir / (f"{rmsk_result_prefix}_{cls}.tsv") for cls in rmsk_classes],
        [
            rmsk_results_dir / (f"{rmsk_result_prefix}_{cls}_{fam}.tsv")
            for cls, fams in rmsk_classes.items()
            for fam in fams
        ],
    conda:
        str(envs_dir / "bedtools.yml")
    log:
        rmsk_results_dir / "rmsk.log",
    params:
        prefix=rmsk_result_prefix,
        filt=lookup_global_chr_filter(config),
    script:
        str(scripts_dir / "get_repeat_masker_features.py")
