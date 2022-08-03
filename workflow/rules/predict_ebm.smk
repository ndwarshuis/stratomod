predict_dir = ebm_dir / "predictions" / "{predict_key}"
predict_log_dir = predict_dir / "log"

rule add_prediction_annotations:
    input:
        variants=rules.concat_tsv_files.output,
        tsvs=[
            expand(
                rules.get_homopolymers.output,
                bases=config["features"]["homopolymers"]["bases"],
                allow_missing=True,
            ),
            rules.get_repeat_masker_classes.output,
            rules.get_tandem_repeats.output,
            rules.get_mappability_high_src.output,
            rules.subtract_high_from_low_mappability.output,
            rules.get_segdups.output,
        ],
    output:
        annotated_input_dir / "{filter_key}.tsv",
    conda:
        str(envs_dir / "bedtools.yml")
    log:
        annotated_input_dir / "{filter_key}.log",
    benchmark:
        annotated_input_dir / "{filter_key}.bench"
    resources:
        mem_mb=attempt_mem_gb(32),
    script:
        str(scripts_dir / "annotate.py")

rule postprocess_predict_x:
    input:
        lambda wildcards: expand(
            rules.add_prediction_annotations.output,
            allow_missing=True,
            input_key=wildcards.input_keys.split(input_delim),
        ),
    output:
        df=predict_dir / "predict_x.tsv",
    params:
        config=lambda wildcards: lookup_ebm_run(wildcards),
    log:
        predict_log_dir / "postprocess.log",
    benchmark:
        predict_log_dir / "postprocess.bench"
    resources:
        mem_mb=attempt_mem_gb(8),
    # TODO somehow need to somehow make this script not freak out when there are
    # no labels
    script:
        str(scripts_dir / "postprocess.py")

rule run_prediction:
    input:
        model = rules.train_ebm.output.model
        predict_x = rules.postprocess_predict_x.output
    output:
        predictions = predict_dir / "predictions.csv"
        explanations = predict_dir / "explanation.csv"
    conda:
        str(envs_dir / "ebm.yml")
    # TODO will likely need this later to get the features for headers
    # params:
    #     config=lambda wildcards: lookup_ebm_run(wildcards),
    log:
        predict_log_dir / "model.log",
    resources:
        mem_mb=attempt_mem_gb(2),
    benchmark:
        predict_log_dir / "model.bench"
    script:
        str(scripts_dir / "annotate.py")
