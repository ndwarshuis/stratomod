import re
from more_itertools import flatten, duplicates_everseen
from itertools import product
from functools import partial

# ------------------------------------------------------------------------------
# too useful...


def fmt_strs(ss):
    return ", ".join(ss)


# ------------------------------------------------------------------------------
# assertions


def assert_empty(xs, msg):
    assert len(xs) == 0, f"{msg}: {fmt_strs(xs)}"


def assert_no_dups(xs, msg):
    assert_empty(set(duplicates_everseen(xs)), msg)


# ------------------------------------------------------------------------------
# Validation
#
# These functions are necessary to validate properties of the config file(s)
# that are not possible using JSON schema validation.
#
# In particular:
# - For each in 'ebm_runs', make sure there are no duplicated feature names
#   (especially considering that features can be renamed on-the-fly with
#   'alt_name' and make sure 'alt_name' features have the same prefix as
#   their parent feature name if they exist
# - For each in 'ebm_runs', make sure explicitly named interaction terms are
#   also in the feature set
# - Make sure input files in for each in 'ebm_runs' are in the 'inputs' section
# - Ensure all keys in input section (the keys under 'inputs' and the keys
#   under each 'test') are unique


def validate_inputs(config):
    train = input_train_keys(config)
    test = input_test_keys(config)
    assert_no_dups(train + test, "Duplicate input keys found")


# TODO make unittests for these
def validate_ebm_features(config):
    def assert_1(run_name, feature_list, feature):
        assert (
            feature in feature_list
        ), f"Interaction {feature} not found in feature set for run {run_name}"

    def assert_N(run_name, feature_list, ints):
        for i in ints:
            assert_1(run_name, feature_list, i)

    def flatten_features(fs):
        return [k if v["alt_name"] is None else v["alt_name"] for k, v in fs.items()]

    flat = [
        (k, ints, v["features"])
        for k, v in config["ebm_runs"].items()
        if "interactions" in v and isinstance(ints := v["interactions"], list)
    ]

    for run_name, ints, features in flat:
        # test that all feature names are valid
        valid_features = set(all_feature_names(config))
        fs = set(list(features))
        assert (
            fs <= valid_features
        ), f"Invalid features: {fmt_strs(fs - valid_features)}"

        for k, v in features.items():
            # test feature alt names
            alt = v["alt_name"]
            if alt is not None:
                prefix = re.match("^[^_]+", k)[0]
                alt_prefix = re.match("^[^_]+", alt)[0]
                assert alt_prefix == prefix, f"Alt prefix must match for {k}"

            # test truncation values
            truncate = v["visualization"]["truncate"]
            tlower = truncate["lower"]
            tupper = truncate["upper"]
            if tlower is not None and tupper is not None:
                assert tlower < tupper, f"Non-positive truncation for {k}"

        # test duplicate alt names
        flat_feature_names = flatten_features(features)
        assert_no_dups(flat_feature_names, "Duplicated features")

        # test for matching interaction terms
        for i in ints:
            check = assert_N if isinstance(i, list) else assert_1
            check(run_name, flat_feature_names, i)


def validate_ebm_inputs(config):
    def assert_keys_exist(what, get_input_keys, get_ebm_keys):
        iks = get_input_keys(config)
        eks = [*flatten(get_ebm_keys(e) for e in config["ebm_runs"].values())]
        assert_empty(
            set(iks) - set(eks),
            f"EBM {what} keys not found in input keys",
        )

    assert_keys_exist("train", input_train_keys, ebm_run_train_keys)
    assert_keys_exist("test", input_test_keys, ebm_run_test_keys)


# ------------------------------------------------------------------------------
# resources


def attempt_mem_gb(mem_gb):
    # double initial memory on each attempt
    return lambda wildcards, attempt: mem_gb * 1000 * 2 ** (attempt - 1)


# ------------------------------------------------------------------------------
# "macros" for rule expansion


def expand_rules(rules, names, attr):
    return [p for r in names for p in getattr(getattr(rules, r), attr)]


# ------------------------------------------------------------------------------
# wildcard sets (for building targets)


def input_set(config, attr):
    return (
        set(
            flatten(
                [t[attr], *nv["train"][attr]]
                for nk, nv in config["inputs"].items()
                for t in nv["test"].values()
            )
        ),
    )


# ------------------------------------------------------------------------------
# global lookup


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


# ------------------------------------------------------------------------------
# input lookup


def lookup_inputs(config):
    return config["inputs"]


def lookup_train(config, train_key):
    return config["inputs"][train_key]


def lookup_all_train(config):
    return [(k, v["train"]) for k, v in lookup_inputs(config).items()]


def lookup_all_test(config):
    return [ts for i in lookup_inputs(config).values() for ts in i["test"].items()]


def lookup_test(config, test_key):
    return dict(lookup_all_test(config))[test_key]


def input_train_keys(config):
    return [*lookup_inputs(config)]


def input_test_keys(config):
    return [i[0] for i in lookup_all_test(config)]


def flat_inputs_(config):
    return flatten(
        [(k, v["train"]), *v["test"].items()] for k, v in config["inputs"].items()
    )


def flat_input_names(config):
    return [i[0] for i in flat_inputs_(config)]


def flat_inputs(config):
    return dict(
        flatten(
            [(k, v["train"]), *v["test"].items()] for k, v in config["inputs"].items()
        )
    )


def lookup_train_test_input(config, input_key):
    return flat_inputs(config)[input_key]


def flat_chr_filters(config):
    return dict(
        flatten(
            [(k, v["chr_filter"]), *[(t, v["chr_filter"]) for t in v["test"]]]
            for k, v in config["inputs"].items()
        )
    )


def lookup_global_chr_filter(config):
    return [*set(flatten(i["chr_filter"] for i in config["inputs"].values()))]


# ------------------------------------------------------------------------------
# ebm run lookup


def lookup_ebm_run(config, run_key):
    return config["ebm_runs"][run_key]


def ebm_run_train_keys(ebm_run):
    return [*flatten([[*i] for i in ebm_run["inputs"]])]


def ebm_run_test_keys(ebm_run):
    return [*flatten([[*ts] for i in ebm_run["inputs"] for ts in i.values()])]


# ------------------------------------------------------------------------------
# feature naming lookup (enter if you dare)

# NOTE: The reason this is so convoluted is because I wanted to have a way to
# test if a given feature set is valid, which means we need to know a priori
# what the total possible feature set is. If this seems unjustified, imagine
# running a gigantic configuration on a cluster, only to have it fail at the
# training step after several hours because a feature was named incorrectly.


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
