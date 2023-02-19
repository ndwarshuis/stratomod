annotations_dir = "annotations"
annotations_src_dir = resources_dir / annotations_dir / cfg.all_wildcards["ref_key"]
annotations_tsv_dir = results_dir / annotations_dir / cfg.all_wildcards["refset_key"]
annotations_log_dir = log_dir / annotations_dir / cfg.all_wildcards["refset_key"]

# All these files from from here:
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

# NOTE: In all cases make no assumptions about sorting or chromosome naming
# standardization. These are all dealt during the initial steps in processing.


include: "genome.smk"
include: "repeat_masker.smk"
include: "homopolymers.smk"
include: "mappability.smk"
include: "tandem_repeats.smk"
include: "segdups.smk"
