from scripts.python.common.config import attempt_mem_gb

tr_dir = "tandem_repeats"


use rule download_labeled_query_vcf as download_tandem_repeats with:
    output:
        config.annotation_resource_dir(tr_dir) / "simple_repeats.txt.gz",
    params:
        src=lambda w: config.references[w.ref_key].annotations.simreps.src,
    localrule: True


# NOTE sorting is done internally by the script
rule get_tandem_repeats:
    input:
        src=partial(
            expand_refkey_from_refsetkey,
            rules.download_tandem_repeats.output,
        ),
        genome=rules.get_genome.output,
    output:
        ensure(
            config.annotation_dir(tr_dir, log=False) / "tandem_repeats.tsv.gz",
            non_empty=True,
        ),
    conda:
        config.env_file("bedtools")
    log:
        config.annotation_dir(tr_dir, log=True) / tr_dir / "tandem_repeats.log",
    benchmark:
        config.annotation_dir(tr_dir, log=True) / "tandem_repeats.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        config.python_script("bedtools/get_tandem_repeat_features.py")
