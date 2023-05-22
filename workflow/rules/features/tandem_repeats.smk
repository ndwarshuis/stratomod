tr_dir = "tandem_repeats"


use rule download_mappability_high as download_tandem_repeats with:
    output:
        config.features_src_dir(log=False) / tr_dir / "simple_repeats.txt.gz",
    log:
        config.features_src_dir(log=True) / tr_dir / "download.log",
    params:
        src=lambda w: config.references[w.ref_key].feature_data.tandem_repeats.src,
    localrule: True


rule get_tandem_repeats:
    input:
        src=partial(
            expand_refkey_from_refsetkey,
            rules.download_tandem_repeats.output,
        ),
        genome=rules.fasta_to_genome.output,
    output:
        ensure(
            config.features_res_dir(log=False) / tr_dir / "tandem_repeats.tsv.gz",
            non_empty=True,
        ),
    conda:
        "../../envs/bio.yml"
    log:
        config.features_res_dir(log=True) / tr_dir / "tandem_repeats.log",
    benchmark:
        config.features_res_dir(log=True) / tr_dir / "tandem_repeats.bench"
    script:
        "../../scripts/python/bio/get_tandem_repeat_features.py"
