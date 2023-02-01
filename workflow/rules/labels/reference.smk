from os.path import dirname
from functools import partial
from scripts.python.common.config import refsetkey_to_chr_filter
from scripts.python.common.functional import compose

ref_resources_dir = resources_dir / "reference" / all_wildcards["ref_key"]
refset_dir = results_dir / all_wildcards["refset_key"]

################################################################################
# download reference

refset_ref_dir = refset_dir / "reference"


rule download_ref_sdf:
    output:
        directory(ref_resources_dir / "sdf"),
    params:
        url=partial(refkey_to_ref_wc, ["sdf"]),
    conda:
        envs_path("utils.yml")
    shell:
        """
        mkdir {output} && \
        curl -Ss {params.url} | \
        bsdtar -xf - \
        --directory {output} \
        --strip-components=1
        """


rule sdf_to_fasta:
    input:
        partial(expand_refkey_from_refsetkey, rules.download_ref_sdf.output),
    output:
        refset_ref_dir / "ref.fa",
    params:
        filt=compose(" ".join, refsetkey_to_chr_filter_wc),
    conda:
        envs_path("rtg.yml")
    benchmark:
        refset_ref_dir / "ref.bench"
    shell:
        """
        rtg sdf2fasta \
        -Z --line-length=70 -n \
        -i {input} \
        -o {output} {params.filt}
        """


################################################################################
# download stratifications

# so far the only use for the stratifications here is to remove MHC since some
# benchmarks don't have VAF/DP here and removing only from the benchmark would
# produce automatic FPs


rule download_mhc_strat:
    output:
        ref_resources_dir / "strats" / "mhc.bed.gz",
    params:
        url=partial(refkey_to_ref_wc, ["strats", "mhc"]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"
