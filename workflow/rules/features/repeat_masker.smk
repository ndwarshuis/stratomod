from scripts.python.common.config import attempt_mem_gb

rmsk_dir = "repeat_masker"

rmsk_res = config.features_res_dir(log=False) / rmsk_dir
rmsk_log = config.features_res_dir(log=True) / rmsk_dir
rmsk_wc_constraint = "[A-Za-z0-9-]+"


def class_file(ext):
    return f"{{rmsk_class}}.{ext}"


def family_file(ext):
    return f"{{rmsk_class}}_{{rmsk_family}}.{ext}"


use rule download_mappability_high as download_repeat_masker with:
    output:
        config.features_src_dir(log=False) / rmsk_dir / "repeat_masker.txt.gz",
    log:
        config.features_src_dir(log=True) / rmsk_dir / "download.log",
    params:
        src=lambda w: config.refkey_to_feature_data(w.ref_key).repeat_masker.src,
    localrule: True


rule get_repeat_masker_classes:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_repeat_masker.output),
    output:
        rmsk_res / class_file("tsv.gz"),
    log:
        rmsk_log / class_file("log"),
    benchmark:
        rmsk_log / class_file("bench")
    conda:
        "../../envs/bio.yml"
    resources:
        mem_mb=attempt_mem_gb(2),
    wildcard_constraints:
        rmsk_class=rmsk_wc_constraint,
    script:
        "../../scripts/python/bio/get_repeat_masker_features.py"


use rule get_repeat_masker_classes as get_repeat_masker_families with:
    output:
        rmsk_res / family_file("tsv.gz"),
    log:
        rmsk_log / family_file("log"),
    benchmark:
        rmsk_log / family_file("bench")
    wildcard_constraints:
        rmsk_class=rmsk_wc_constraint,
        family_class=rmsk_wc_constraint,


def rmsk_targets(ref_key):
    cfs = config.references[ref_key].feature_data.repeat_masker.class_families
    cs, fs = zip(*[(c, f) for c, fs in cfs.items() for f in fs])
    return expand(
        rules.get_repeat_masker_classes.output,
        allow_missing=True,
        rmsk_class=list(cfs),
    ) + expand(
        rules.get_repeat_masker_families.output,
        zip,
        allow_missing=True,
        rmsk_class=cs,
        rmsk_family=fs,
    )
