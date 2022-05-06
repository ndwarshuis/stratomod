from more_itertools import flatten
from itertools import product
from functools import partial


def lookup_config(config, *keys):
    k = keys[0]
    ks = keys[1:]
    return config[k] if len(ks) == 0 else lookup_config(config[k], *ks)


def lookup_annotations(config):
    # TODO sometime in the future if/when we use more than just GRCh38, don't
    # hard-code this
    return lookup_config(
        config,
        "resources",
        "references",
        "GRCh38",
        "annotations",
    )


def lookup_global_chr_filter(config):
    return [*set(flatten(i["chr_filter"] for i in config["inputs"].values()))]


def lookup_bed_cols(config):
    return config["features"]["index"]


def bed_cols_ordered(bed_cols):
    return [bed_cols["chr"], bed_cols["start"], bed_cols["end"]]


def bed_cols_indexed(indices, bed_cols):
    return dict(zip(indices, bed_cols_ordered(bed_cols)))


def lookup_bed_cols_ordered(config):
    return bed_cols_ordered(lookup_bed_cols(config))


def fmt_feature(prefix, rest):
    return f"{prefix}_{rest}"


def fmt_vcf_feature(config, which):
    fconf = config["features"]["vcf"]
    return fmt_feature(fconf["prefix"], fconf["columns"][which])


def vcf_feature_names(config):
    return [
        *map(
            lambda f: fmt_vcf_feature(config, f),
            config["features"]["vcf"]["columns"],
        )
    ]


def fmt_mappability_feature(config, which):
    fconf = config["features"]["mappability"]
    return fmt_feature(fconf["prefix"], fconf["suffixes"][which])


def mappability_feature_names(config):
    return [
        *map(
            lambda f: fmt_mappability_feature(config, f),
            config["features"]["mappability"]["suffixes"],
        ),
    ]


def fmt_homopolymer_feature(config, bases, which):
    fconf = config["features"]["homopolymers"]
    return fmt_feature(fconf["prefix"], f"{bases}_{fconf['suffixes'][which]}")


def homopolymer_feature_names(config):
    fconf = config["features"]["homopolymers"]
    return [
        fmt_homopolymer_feature(config, b, s)
        for (b, s) in product(fconf["bases"], fconf["suffixes"])
    ]


def fmt_repeat_masker_feature(config, grp, fam=None):
    fconf = config["features"]["repeat_masker"]
    rest = grp if fam is None else f"{grp}_{fam}"
    # TODO this is hardcoded for now
    suffix = fconf["suffixes"]["len"]
    return fmt_feature(fconf["prefix"], f"{rest}_{suffix}")


def repeat_masker_feature_names(config):
    fconf = config["features"]["repeat_masker"]
    fmt = partial(fmt_repeat_masker_feature, config)
    return [
        *[fmt(c) for c in fconf["classes"]],
        *[fmt(c, f) for c, fs in fconf["classes"].items() for f in fs],
    ]


def fmt_count_feature(prefix):
    return fmt_feature(prefix, "count")


def fmt_merged_feature(prefix, middle, op):
    return fmt_feature(prefix, f"{middle}_{op}")


def merged_feature_names(prefix, names, ops):
    return [
        *[fmt_merged_feature(prefix, n, o) for n, o in product(names, ops)],
        fmt_count_feature(prefix),
    ]


def segdup_feature_names(config):
    fconf = config["features"]["segdups"]
    prefix = fconf["prefix"]
    return merged_feature_names(
        prefix,
        fconf["columns"].values(),
        fconf["operations"],
    )


def fmt_tandem_repeat_base(config, bases):
    bs_prefix = config["features"]["tandem_repeats"]["bases_prefix"]
    return f"{bs_prefix}_{bases}"


def tandem_repeat_feature_names(config):
    fconf = config["features"]["tandem_repeats"]
    prefix = fconf["prefix"]
    # TODO weirdly hardcoded in several places
    bases = ["A", "T", "G", "C", "AT", "GC"]
    bs = [fmt_tandem_repeat_base(config, b) for b in bases]
    cs = fconf["columns"].values()
    return [
        *merged_feature_names(prefix, [*cs, *bs], fconf["operations"]),
        *[fmt_feature(prefix, o) for o in fconf["other"].values()],
    ]


def all_feature_names(config):
    return [
        *vcf_feature_names(config),
        *mappability_feature_names(config),
        *homopolymer_feature_names(config),
        *repeat_masker_feature_names(config),
        *segdup_feature_names(config),
        *tandem_repeat_feature_names(config),
    ]
