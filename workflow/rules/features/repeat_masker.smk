from scripts.python.common.config import attempt_mem_gb

rmsk_dir = "repeat_masker"
rmsk_classes = config.feature_names.repeat_masker.classes

rmsk_res = config.features_res_dir(log=False) / rmsk_dir
rmsk_log = config.features_res_dir(log=True) / rmsk_dir


use rule download_mappability_high as download_repeat_masker with:
    output:
        config.features_src_dir(log=False) / rmsk_dir / "repeat_masker.txt.gz",
    log:
        config.features_src_dir(log=True) / rmsk_dir / "download.log",
    params:
        src=lambda w: config.refkey_to_feature_data(w.ref_key).repeat_masker.src,
    localrule: True


# use keyed outputs to allow easier parsing in the script
rule get_repeat_masker_classes:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_repeat_masker.output),
    output:
        **{
            cls: ensure(rmsk_res / (f"{cls}.tsv.gz"), non_empty=True)
            for cls in rmsk_classes
        },
        **{
            f"{cls}_{fam}": ensure(rmsk_res / (f"{cls}_{fam}.tsv.gz"), non_empty=True)
            for cls, fams in rmsk_classes.items()
            for fam in fams
        },
    conda:
        "../../envs/bio.yml"
    log:
        rmsk_log / "rmsk.log",
    benchmark:
        rmsk_log / "rmsk.bench"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        "../../scripts/python/bio/get_repeat_masker_features.py"
