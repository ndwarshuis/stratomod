from scripts.common.config import lookup_config

inputs_dir = resources_dir / "inputs"

label_dir = results_dir / "labels" / "{input_key}"
rtg_dir = label_dir / "rtg"

labels = ["fp", "fn", "tp"]


include: "download_resources.smk"


################################################################################
# VCF preprocessing


# TODO this is (probably) just for DV VCFs
rule preprocess_vcf:
    input:
        lambda wildcards: inputs_dir
        / lookup_config(config, "inputs", wildcards.input_key, "vcf"),
    output:
        label_dir / "query.vcf.gz",
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        """
        gunzip -c {input} | \
        sed -e '/.RefCall./ s/\.\/\./0\/1/g' | \
        sed -e '/.RefCall./ s/0\/0/0\/1/g' | \
        bgzip -c > {output}
        """


rule index_vcf:
    input:
        rules.preprocess_vcf.output,
    output:
        label_dir / "query.vcf.gz.tbi",
    conda:
        str(envs_dir / "samtools.yml")
    shell:
        "tabix -p vcf {input}"


################################################################################
# VCF -> tsv


# rtg won't output to a directory that already exists, so do this weird temp
# file thing
# TODO add option to switch of the "--ref-overlap --all-records" thingy
# TODO --all-records = use all records including those that fail, make an option
# for this


def lookup_input(wildcards, *args):
    return lookup_config(config, "inputs", wildcards.input_key, *args)


def get_truth_inputs(wildcards):
    out = rules.get_bench.output
    return {
        key: expand(path, bench_key=lookup_input(wildcards, "benchmark"))
        for key, path in [
            ("truth_vcf", out.vcf),
            ("truth_bed", out.bed),
            ("truth_tbi", out.tbi),
        ]
    }


rule get_vcf_labels:
    input:
        unpack(get_truth_inputs),
        query_vcf=rules.preprocess_vcf.output,
        sdf=lambda wildcards: expand(
            rules.get_ref_sdf.output,
            ref_key=lookup_input(wildcards, "ref"),
        ),
        # not used on CLI but still needed
        query_tbi=rules.index_vcf.output,
    output:
        [rtg_dir / ("%s.vcf.gz" % lbl) for lbl in labels]
    conda:
        str(envs_dir / "rtg.yml")
    params:
        extra="--ref-overlap --all-records",
        tmp_dir="/tmp/vcfeval",
        output_dir=lambda _, output: Path(output[0]).parent,
    shell:
        """
        rtg vcfeval {params.extra} \
            -b {input.truth_vcf} \
            -e {input.truth_bed} \
            -c {input.query_vcf} \
            -o {params.tmp_dir} \
        -t {input.sdf}

        mv {params.tmp_dir}/* {params.output_dir}

        rm -r {params.tmp_dir}
        """


rule unzip_vcf_labels:
    input:
        rtg_dir / "{label}.vcf.gz",
    output:
        label_dir / "{label}.vcf",
    shell:
        "gunzip {input} -c > {output}"


rule parse_label_vcf:
    input:
        rules.unzip_vcf_labels.output,
    output:
        label_dir / "{filter_key}_{label}.tsv",
    shell:
        """
        python \
        workflow/scripts/parse_vcf_to_bed_ebm.py \
        --type {wildcards.filter_key} \
        --label {wildcards.label} \
        --input {input} \
        --output {output}
        """

rule concat_tsv_files:
    input:
        expand(rules.parse_label_vcf.output, label=labels, allow_missing=True)
    output:
        label_dir / "{filter_key}_labeled.tsv",
    run:
        import pandas as pd
        import workflow.scripts.common.tsv
        import workflow.scripts.common.bed
        # use pandas here since it will more reliably account for headers
        df = pd.concat([read_tsv(i, header=0) for i in input])
        write_tsv(output, sort_bed_numerically(df))
    # shell:
    #     """
    #     tail -n+2 {input[0]} | \
    #     cat {input[1]} - | \
    #     cat {input[2]} - | \
    #     python workflow/scripts/sort_and_filter_bed.py --header \
    #     > {output}
    #     """


## TODO add filtering rules here if we wish
