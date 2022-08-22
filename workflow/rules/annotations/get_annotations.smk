annotations_src_dir = resources_dir / "annotations"
annotations_tsv_dir = results_dir / "annotations"

# All these files from from here:
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/


include: "genome.smk"
include: "repeat_masker.smk"
include: "homopolymers.smk"
include: "mappability.smk"
include: "tandem_repeats.smk"
include: "segdups.smk"
