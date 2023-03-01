import scripts.python.common.config as cfg


wildcard_constraints:
    **cfg._constraints,


################################################################################
# wildcard formating


def wildcard_ext(key, ext):
    return f"{cfg.all_wildcards[key]}.{ext}"


def wildcard_format(format_str, *keys):
    return format_str.format(*[cfg.all_wildcards[k] for k in keys])


def wildcard_format_ext(format_str, keys, ext):
    return wildcard_format(f"{format_str}.{ext}", *keys)


################################################################################
# config lookup


def refsetkey_to_refkey_wc(wildcards):
    return config.refsetkey_to_refkey(wildcards.refset_key)


def refkey_to_ref_wc(args, wildcards):
    return config.refkey_to_ref(wildcards.ref_key)


def inputkey_to_input_wc(args, wildcards):
    return config.inputkey_to_input(wildcards.input_key)


def inputkey_to_chr_filter_wc(wildcards):
    filt = config.inputkey_to_chr_filter(wildcards.input_key)
    return "\\|".join(f"{f}\\b" for f in filt)


def refsetkey_to_chr_indices_wc(wildcards):
    return config.refsetkey_to_chr_indices(wildcards.refset_key)


# def refsetkey_to_benchmark_wc(key, wildcards):
#     return config.refsetkey_to_benchmark([key], wildcards.bench_key, wildcards.refset_key)


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
