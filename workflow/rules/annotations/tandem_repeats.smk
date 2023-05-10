from scripts.python.common.config import attempt_mem_gb

tr_dir = "tandem_repeats"


use rule download_mappability_high as download_tandem_repeats with:
    output:
        config.annotation_src_dir(log=False) / tr_dir / "simple_repeats.txt.gz",
    log:
        config.annotation_src_dir(log=True) / tr_dir / "download.log",
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
        genome=rules.fasta_to_genome.output,
    output:
        ensure(
            config.annotation_res_dir(log=False) / tr_dir / "tandem_repeats.tsv.gz",
            non_empty=True,
        ),
    conda:
        "../../envs/bio.yml"
    log:
        config.annotation_res_dir(log=True) / tr_dir / "tandem_repeats.log",
    benchmark:
        config.annotation_res_dir(log=True) / tr_dir / "tandem_repeats.bench"
    resources:
        mem_mb=attempt_mem_gb(1),
    script:
        "../../scripts/python/bio/get_tandem_repeat_features.py"
