import re
from typing import Any
from more_itertools import flatten, duplicates_everseen, unzip
from itertools import product, filterfalse
from functools import partial
from .functional import maybe, compose

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
            set(eks) - set(iks),
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
# wildcard sets (for building targets)


def all_benchkeys(config: dict) -> dict[str, list[str]]:
    bks = [
        (inputkey_to_refkey(config, k), b)
        for k, v in flat_inputs(config).items()
        if (b := v["benchmark"]) is not None
    ]
    return {k: list(v) for k, v in zip(["bench_key", "ref_key"], unzip(bks))}


def input_set(config, key):
    return set(
        filterfalse(lambda x: not x, [x[key] for x in flat_inputs(config).values()])
    )
    # flatten(
    #     [
    #         *[a for t in nv["test"].values() if (a := t[attr]) is not None],
    #         nv["train"][attr],
    #     ]
    #     for nk, nv in config["inputs"].items()
    # )
    # )


# ------------------------------------------------------------------------------
# global lookup


def walk_dict(d, keys):
    k = keys[0]
    ks = keys[1:]
    dnext = d[k]
    return dnext if len(ks) == 0 else walk_dict(dnext, ks)


def lookup_config(config, keys):
    return walk_dict(config, keys)


def chr_index_to_str(i: int) -> str:
    return "X" if i == 23 else ("Y" if i == 24 else str(i))


def chr_indices_to_name(prefix: str, xs: list[int]):
    return [f"{prefix}{chr_index_to_str(i)}" for i in xs]


# ------------------------------------------------------------------------------
# lookup from ref_key


def refkey_to_ref(config: dict, args: list[str], ref_key: str) -> Any:
    return lookup_config(config, ["references", ref_key, *args])


def refkey_to_benchmark(
    config: dict,
    args: list[str],
    bench_key: str,
    ref_key: str,
) -> Any:
    return refkey_to_ref(config, ["benchmarks", bench_key, *args], ref_key)


# ------------------------------------------------------------------------------
# lookup from refset_key


def refsetkey_to_benchmark(
    config: dict,
    args: list[str],
    bench_key: str,
    refset_key: str,
) -> Any:
    return compose(
        partial(refkey_to_benchmark, config, args, bench_key),
        partial(refsetkey_to_refkey, config),
    )(refset_key)


def refsetkey_to_refset(config: dict, args: list[str], refset_key: str):
    return lookup_config(config, ["reference_sets", refset_key, *args])


def refsetkey_to_refkey(config, refset_key):
    return refsetkey_to_refset(config, ["ref"], refset_key)


def refsetkey_to_ref(config: dict, args: list[str], refset_key: str) -> Any:
    return compose(
        partial(refkey_to_ref, config, args),
        partial(inputkey_to_refkey, config, refset_key),
    )


def refsetkey_to_chr_indices(config: dict, refset_key: str) -> list[int]:
    f = refsetkey_to_refset(config, ["chr_filter"], refset_key)
    return list(range(1, 25)) if len(f) == 0 else f


def refsetkey_to_chr_prefix(config: dict, refset_key: str) -> str:
    return compose(
        partial(refkey_to_ref, config, ["chr_prefix"]),
        partial(refsetkey_to_refkey, config),
    )(refset_key)


def refsetkey_to_chr_filter(config: dict, refset_key: str) -> list[str]:
    indices = refsetkey_to_chr_indices(config, refset_key)
    prefix = refsetkey_to_ref(config, ["chr_prefix"], refset_key)
    return chr_indices_to_name(prefix, indices)


# ------------------------------------------------------------------------------
# lookup from input_key

# This is tricky because 99% of the time I want to think of train and test
# inputs as being part of a flat list; therefore, "input_key" means "train_key
# or test_key" (most of the time)


def inputkey_to_input(config: dict, args: list[str], input_key: str):
    return walk_dict(flat_inputs(config)[input_key], args)


def inputkey_to_shared(config: dict, input_key: str, shared_key: str) -> Any:
    """
    Given an input key, return a value indicated by 'shared_key', which is a
    key common to all input dictionaries.
    """
    return dict(
        flatten(
            [
                (train_key, v[shared_key]),
                *[(test_key, v[shared_key]) for test_key in v["test"]],
            ]
            for train_key, v in config["inputs"].items()
        )
    )[input_key]


def inputkey_to_refkey(config: dict, input_key: str):
    return compose(
        partial(refsetkey_to_refkey, config),
        partial(inputkey_to_refsetkey, config),
    )(input_key)


def inputkey_to_refsetkey(config: dict, input_key: str):
    return inputkey_to_shared(config, input_key, "refset")


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


def all_refsetkeys(config):
    return set(v["refset"] for v in config["inputs"].values())


def all_refkeys(config):
    return set(map(partial(refsetkey_to_refkey, config), all_refsetkeys(config)))


def flat_input_names(config):
    return [i[0] for i in flat_inputs_(config)]


# TODO this isn't a well defined function in terms of types since the test
# dictionary has several optional fields and the train does not.
def flat_inputs(config: dict) -> dict[str, Any]:
    """Return a dictionary of all inputs in the config."""
    return dict(
        flatten(
            [(k, v["train"]), *v["test"].items()] for k, v in config["inputs"].items()
        )
    )


def inputkey_to_chr_prefix(config: dict, input_key: str) -> str:
    input_prefix = inputkey_to_input(config, ["chr_prefix"], input_key)
    if input_prefix is None:
        return compose(
            partial(refsetkey_to_chr_prefix, config),
            partial(inputkey_to_refsetkey, config),
        )(input_key)
    else:
        return input_prefix


def inputkey_to_chr_filter(config: dict, input_key: str) -> list[str]:
    prefix = inputkey_to_chr_prefix(config, input_key)
    return compose(
        partial(chr_indices_to_name, prefix),
        partial(refsetkey_to_chr_indices, config),
        partial(inputkey_to_refsetkey, config),
    )(input_key)


def inputkey_to_bench_correction(config: dict, key: str, input_key: str) -> bool:
    bench_key = inputkey_to_input(config, ["benchmark"], input_key)
    return compose(
        partial(refkey_to_benchmark, config, ["corrections", key], bench_key),
        partial(refsetkey_to_refkey, config),
        partial(inputkey_to_refsetkey, config),
    )(input_key)


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


def lookup_raw_index(config):
    return config["features"]["raw_index"]


def lookup_bed_cols(config):
    return config["features"]["bed_index"]


def bed_cols_ordered(bed_cols):
    return [bed_cols["chr"], bed_cols["start"], bed_cols["end"]]


def bed_cols_indexed(indices, bed_cols):
    return dict(zip(indices, bed_cols_ordered(bed_cols)))


def lookup_bed_cols_ordered(config):
    return bed_cols_ordered(lookup_bed_cols(config))


def lookup_all_index_cols(config):
    return [lookup_raw_index(config), *bed_cols_ordered(lookup_bed_cols(config))]


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
    rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
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
