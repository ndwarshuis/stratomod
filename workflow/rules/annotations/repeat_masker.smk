from scripts.python.common.config import attempt_mem_gb

rmsk_dir = "repeat_masker"
rmsk_classes = config.feature_names.repeat_masker.classes

rmsk_res = config.annotation_dir(rmsk_dir, log=False)
rmsk_log = config.annotation_dir(rmsk_dir, log=True)


rule download_repeat_masker:
    output:
        config.annotation_resource_dir(rmsk_dir) / "repeat_masker.txt.gz",
    params:
        url=lambda wildcards: config.refkey_to_annotations(
            wildcards.ref_key
        ).repeat_masker.url,
    conda:
        config.env_file("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


# use keyed outputs to allow easier parsing in the script
rule get_repeat_masker_classes:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_repeat_masker.output),
    output:
        **{
            cls: ensure(rmsk_res / (f"{cls}.tsv.gz"), non_empty=True)
            for cls in rmsk_classes.dict()
        },
        **{
            f"{cls}_{fam}": ensure(rmsk_res / (f"{cls}_{fam}.tsv.gz"), non_empty=True)
            for cls, fams in rmsk_classes.dict().items()
            for fam in fams
        },
    conda:
        config.env_file("bedtools")
    log:
        rmsk_log / "rmsk.log",
    benchmark:
        rmsk_log / "rmsk.bench"
    resources:
        mem_mb=attempt_mem_gb(2),
    script:
        config.python_script("bedtools/get_repeat_masker_features.py")
