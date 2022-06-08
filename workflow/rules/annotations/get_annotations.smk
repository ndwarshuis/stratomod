annotations_src_dir = resources_dir / "annotations"
annotations_tsv_dir = results_dir / "annotations"

# All these files from from here:
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

# Download genome file (necessary from pretty much all the specific annotation
# rules below). Use the first two columns of this table (chrom and length)
rule get_genome:
    output:
        annotations_src_dir / "genome.txt",
    params:
        url=config["resources"]["references"]["GRCh38"]["genome"]
    conda:
        str(envs_dir / "utils.yml")
    shell:
        """
        curl -Ss {params.url} | \
        gunzip -c | \
        cut -f1,2 | \
        sed -n '/^chr\([0-9XY][[:space:]]\|[0-9]\{{2\}}[[:space:]]\)/p' | \
        sed 's/^chr//' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' \
        > {output}
        """

include: "repeat_masker.smk"
include: "homopolymers.smk"
include: "mappability.smk"
include: "tandem_repeats.smk"
include: "segdups.smk"
