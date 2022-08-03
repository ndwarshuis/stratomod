predict_dir = ebm_dir / "predictions" / "{predict_key}"
predict_log_dir = predict_dir / "log"

# rule postprocess_predict_x:
#     input:

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
