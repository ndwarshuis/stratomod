import scripts.python.common.config as cfg

# wildcards shall not be used for multiple nodes in a filepath
# TODO there are lots of characters that should be avoided to prevent ambiguity

ALL_LABELS = ["fp", "fn", "tp"]

ALL_FILTERS = ["SNP", "INDEL"]

all_bases = config["features"]["homopolymers"]["bases"]


def alternate_constraint(xs):
    return f"({'|'.join(xs)})"


_constraints = {
    # corresponds to a genome reference
    "ref_key": "[^/]+",
    # corresponds to a reference set (reference + chromosome filter + etc)
    "refset_key": "[^/]+",
    # corresponds to a training vcf
    "train_key": "[^/]+",
    # corresponds to a testing vcf
    "test_key": "[^/]+",
    # refers to either a train_key or test_key (for instances where they are to
    # be treated the same)
    "input_key": "[^/]+",
    # refers to a collection of input data input to an ebm model configuration
    # (composed of multiple train/test keys + associated data)
    "run_key": "[^/]+",
    # refers to a benchmark vcf (within the context of a given reference)
    "bench_key": "[^/]+",
    # refers to an EBM model and its parameters
    "model_key": "[^/]+",
    # refers to the variant type (SNP or INDEL, for now)
    # TODO ...why "filter"? (I can't think of anything better)
    "filter_key": alternate_constraint(ALL_FILTERS),
    # refers to a variant benchmarking label (tp, fp, etc)
    "label": alternate_constraint(ALL_LABELS),
    # refers to a nucleotide base
    "base": alternate_constraint(all_bases),
}


wildcard_constraints:
    **_constraints,


all_wildcards = {k: f"{{{k},{v}}}" for k, v in _constraints.items()}


################################################################################
# wildcard formating


def wildcard_ext(key, ext):
    return f"{all_wildcards[key]}.{ext}"


def wildcard_format(format_str, *keys):
    return format_str.format(*[all_wildcards[k] for k in keys])


def wildcard_format_ext(format_str, keys, ext):
    return wildcard_format(f"{format_str}.{ext}", *keys)


################################################################################
# config lookup


def refsetkey_to_refkey_wc(wildcards):
    return cfg.refsetkey_to_refkey(config, wildcards.refset_key)


def refkey_to_ref_wc(args, wildcards):
    return cfg.refkey_to_ref(config, args, wildcards.ref_key)


def inputkey_to_input_wc(args, wildcards):
    return cfg.inputkey_to_input(config, args, wildcards.input_key)


def inputkey_to_chr_filter_wc(wildcards):
    filt = cfg.inputkey_to_chr_filter(config, wildcards.input_key)
    return "\\|".join(f"{f}\\b" for f in filt)


def refsetkey_to_chr_indices_wc(wildcards):
    return cfg.refsetkey_to_chr_indices(config, wildcards.refset_key)


def refsetkey_to_benchmark_wc(key, wildcards):
    return cfg.refsetkey_to_benchmark(
        config, [key], wildcards.bench_key, wildcards.refset_key
    )


def refkey_to_benchmark_wc(key, wildcards):
    return cfg.refkey_to_benchmark(
        config,
        [key],
        wildcards.bench_key,
        wildcards.ref_key,
    )


################################################################################
# path expansion


def expand_refkey_from_refsetkey(path, wildcards):
    return expand(
        path,
        allow_missing=True,
        ref_key=refsetkey_to_refkey_wc(wildcards),
    )


def expand_refsetkey_from_inputkey(path, wildcards):
    return expand(
        path,
        allow_missing=True,
        refset_key=inputkey_to_refsetkey(config, wildcards.input_key),
    )


def expand_benchkey_from_inputkey(path, wildcards):
    return expand(
        path,
        allow_missing=True,
        bench_key=inputkey_to_input_wc(["benchmark"], wildcards),
    )


def expand_rules(names, attr):
    return [p for r in names for p in getattr(getattr(rules, r), attr)]


def rule_output_suffix(rulename, suffix):
    return f"{getattr(rules, rulename).output[0]}.{suffix}"
