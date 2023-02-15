import re
from typing import (
    Any,
    Sequence,
    Dict,
    List,
    Tuple,
    Union,
    Set,
    Iterable,
    Collection,
    Callable,
    Optional,
    NamedTuple,
    NewType,
)
from more_itertools import flatten, duplicates_everseen, unzip, partition
from itertools import product
from functools import partial
from .functional import maybe, compose
from snakemake.io import expand, InputFiles  # type: ignore
from pydantic import BaseModel
from enum import Enum

RefKey = NewType("RefKey", str)
RefsetKey = NewType("RefsetKey", str)
TestKey = NewType("TestKey", str)
InputKey = NewType("InputKey", str)
BenchKey = NewType("BenchKey", str)
ModelKey = NewType("ModelKey", str)
RunKey = NewType("RunKey", str)
FeatureKey = NewType("FeatureKey", str)

VarKey = NewType("VarKey", str)


class Base(Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"


class FilterKey(Enum):
    SNV = "SNV"
    INDEL = "INDEL"


class Label(Enum):
    TP = "tp"
    FP = "fp"
    FN = "fn"


class PlotType(Enum):
    BAR = "bar"
    STEP = "step"


class FeatureType(Enum):
    CONTINUOUS = "continuous"
    CATEGORICAL = "categorical"


class Binning(Enum):
    UNIFORM = "uniform"
    QUANTILE = "quantile"
    QUANTILE_HUMANIZED = "quantile_humanized"


class BedMergeOp(Enum):
    MIN = "min"
    MAX = "max"
    MEAN = "mean"
    MEDIAN = "median"


class Transform(Enum):
    LOG = "log"
    BINARY = "binary"


class Paths(BaseModel):
    resources: str
    results: str


class RefSet(BaseModel):
    ref: RefKey
    chr_filter: Set[int]


class CatVar(BaseModel):
    levels: Set[str]


class ContVar(BaseModel):
    lower: Optional[float]
    upper: Optional[float]


class FormatFields(BaseModel):
    vaf: Optional[str]
    dp: Optional[str]
    gt: Optional[str]
    gq: Optional[str]


class VCFInput(BaseModel):
    refset: RefsetKey
    chr_prefix: str
    url: Optional[str]
    benchmark: Optional[BenchKey]
    variables: Dict[VarKey, str]
    format_fields: Optional[FormatFields]
    # TODO min = 1
    max_ref: int = 50
    max_alt: int = 50


class TestDataInput(BaseModel):
    input_key: InputKey
    variables: Dict[VarKey, str]


class ModelInput(BaseModel):
    train: Set[InputKey]
    test: Dict[TestKey, TestDataInput]


class EBMMiscParams(BaseModel):
    # TODO this is a fraction/ needs constraint
    downsample: Optional[float] = None


class EBMSplitParams(BaseModel):
    # TODO this is a fraction/ needs constraint
    test_size: float = 0.2
    random_state: Optional[int] = None


# TODO these all need constraints
class EBMClassifierParams(BaseModel):
    max_bins: int = 256
    max_interaction_bins: int = 32
    binning: Binning = Binning.QUANTILE
    outer_bags: int = 8
    inner_bags: int = 0
    learning_rate: float = 0.01
    validation_size: float = 0.15
    early_stopping_rounds: int = 50
    early_stopping_tolerance: float = 0.0001
    max_rounds: int = 5000
    min_samples_leaf: int = 2
    max_leaves: int = 3
    random_state: Optional[int] = None


class EBMSettings(BaseModel):
    misc_parameters: EBMMiscParams = EBMMiscParams()
    split_parameters: EBMSplitParams = EBMSplitParams()
    classifier_parameters: EBMClassifierParams = EBMClassifierParams()


class Truncation(BaseModel):
    lower: Optional[float]
    upper: Optional[float]


class Visualization(BaseModel):
    truncate: Truncation = Truncation()
    plot_type: PlotType = PlotType.STEP
    split_missing: Optional[float] = None


class Feature(BaseModel):
    feature_type: FeatureType = FeatureType.CONTINUOUS
    fill_na: Optional[float] = 0.0
    alt_name: Optional[str] = None
    visualization: Visualization = Visualization()
    transform: Optional[Transform] = None


class EBMRun(BaseModel):
    inputs: Dict[RunKey, ModelInput]
    filter: Set[FilterKey]
    ebm_settings: EBMSettings
    error_labels: Set[Label]
    filtered_are_candidates: bool
    interactions: Union[int, List[str], List[List[str]]]
    features: Dict[FeatureKey, Feature]


class BedFile(BaseModel):
    url: Optional[str]
    chr_prefix: str


class Strats(BaseModel):
    mhc: BedFile


class BenchmarkCorrections(BaseModel):
    strip_IPS: bool


class Benchmark(BaseModel):
    vcf_url: Optional[str]
    bed_url: Optional[str]
    chr_prefix: str
    corrections: BenchmarkCorrections


class Mappability(BaseModel):
    low: BedFile
    high: BedFile


class Annotations(BaseModel):
    mappability: Mappability
    superdups: BedFile
    simreps: BedFile
    repeat_masker: BedFile


class Reference(BaseModel):
    sdf: BedFile
    genome: BedFile
    strats: Strats
    annotations: Annotations
    benchmarks: Dict[str, Benchmark]


class BedIndex(BaseModel):
    chr: str
    start: str
    end: str


class VCFColumns(BaseModel):
    input: str
    qual: str
    filter: str
    info: str
    gt: str
    gq: str
    dp: str
    vaf: str
    len: str


class VCFMeta(BaseModel):
    prefix: str
    columns: VCFColumns


class MapSuffixes(BaseModel):
    low: str
    high: str


class MapMeta(BaseModel):
    prefix: str
    suffixes: MapSuffixes


class HomopolySuffixes(BaseModel):
    len: str
    imp_frac: str


class HomopolyMeta(BaseModel):
    prefix: str
    bases: List[Base]
    suffixes: HomopolySuffixes


class RMSKSuffixes(BaseModel):
    len: str


class RMSKClasses(BaseModel):
    SINE = List[str]
    LINE = List[str]
    LTR = List[str]
    Satellite = List[str]


class RMSKMeta(BaseModel):
    prefix: str
    suffixes: RMSKSuffixes
    classes: RMSKClasses


class SegDupsColumns(BaseModel):
    alignL: str
    fracMatchIndel: str


class SegDupsMeta(BaseModel):
    prefix: str
    columns: SegDupsColumns
    operations: Set[BedMergeOp]


class TandemRepeatColumns(BaseModel):
    period: str
    copyNum: str
    perMatch: str
    perIndel: str
    score: str


class TandemRepeatOther(BaseModel):
    len: str


class TandemRepeatMeta(BaseModel):
    prefix: str
    bases_prefix: str
    operations: Set[BedMergeOp]
    columns: TandemRepeatColumns
    other: TandemRepeatOther


class FeatureMeta(BaseModel):
    label: str
    raw_index: str
    bed_index: BedIndex
    vcf: VCFMeta
    mappability: MapMeta
    homopolymers: HomopolyMeta
    repeat_masker: RMSKMeta
    segdups: SegDupsMeta
    tandem_repeats: TandemRepeatMeta


class StratoMod(BaseModel):
    reference_sets: Dict[RefsetKey, RefSet]
    variables: Dict[VarKey, Union[ContVar, CatVar]]
    inputs: Dict[InputKey, VCFInput]
    ebm_runs: Dict[ModelKey, EBMRun]
    paths: Paths
    references: Dict[RefKey, Reference]
    feature_meta: FeatureMeta


# ------------------------------------------------------------------------------
# too useful...


def fmt_strs(ss: Iterable[str]) -> str:
    return ", ".join(ss)


# ------------------------------------------------------------------------------
# assertions


def assert_empty(xs: Collection, msg: str) -> None:
    assert len(xs) == 0, f"{msg}: {fmt_strs(xs)}"


def assert_no_dups(xs: Collection, msg: str) -> None:
    assert_empty(set(duplicates_everseen(xs)), msg)


def assert_match(pat: str, s: str) -> str:
    res = re.match(pat, s)
    assert res is not None, f"match failed for pattern '{pat}' and query '{s}'"
    return res[0]


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


def validate_inputs(config: StratoMod) -> None:
    inputs = config.inputs
    # train = input_train_keys(config)
    # test = input_test_keys(config)
    assert_no_dups(inputs, "Duplicate input keys found")


# TODO make unittests for these
def validate_ebm_features(config: StratoMod) -> None:
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
        (k, ints, v.features)
        for k, v in config.ebm_runs.items()
        if isinstance(ints := v.interactions, list)
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
            alt = v.alt_name
            if alt is not None:
                prefix = assert_match("^[^_]+", k)
                alt_prefix = assert_match("^[^_]+", alt)
                assert alt_prefix == prefix, f"Alt prefix must match for {k}"

            # test truncation values
            truncate = v.visualization.truncate
            tlower = truncate.lower
            tupper = truncate.upper
            if tlower is not None and tupper is not None:
                assert tlower < tupper, f"Non-positive truncation for {k}"

        # test duplicate alt names
        flat_feature_names = flatten_features(features)
        assert_no_dups(flat_feature_names, "Duplicated features")

        # test for matching interaction terms
        for i in ints:
            check = assert_N if isinstance(i, list) else assert_1
            check(run_name, flat_feature_names, i)


# TODO fixme
# def validate_ebm_inputs(config: StratoMod) -> None:
#     def assert_keys_exist(what, get_input_keys, get_ebm_keys):
#         iks = get_input_keys(config)
#         eks = [*flatten(get_ebm_keys(e) for e in config["ebm_runs"].values())]
#         assert_empty(
#             set(eks) - set(iks),
#             f"EBM {what} keys not found in input keys",
#         )

#     assert_keys_exist("train", input_train_keys, ebm_run_train_keys)
#     assert_keys_exist("test", input_test_keys, ebm_run_test_keys)


# ------------------------------------------------------------------------------
# resources


def attempt_mem_gb(mem_gb: int) -> Callable[[dict, int], int]:
    # double initial memory on each attempt
    return lambda wildcards, attempt: mem_gb * 1000 * 2 ** (attempt - 1)


# ------------------------------------------------------------------------------
# expanding targets


def all_benchkeys(config: StratoMod, target: InputFiles) -> InputFiles:
    rs, bs = unzip(
        (inputkey_to_refkey(config, k), b)
        for k, v in config.inputs.items()
        if (b := v.benchmark) is not None
    )
    return expand(target, zip, ref_key=rs, bench_key=bs)


class RunKeysTrain(NamedTuple):
    model_key: ModelKey
    filter_key: FilterKey
    data_key: RunKey
    input_key: InputKey


class RunKeysTest(NamedTuple):
    model_key: ModelKey
    filter_key: FilterKey
    data_key: RunKey
    test_key: TestKey
    input_key: InputKey


def lookup_run_sets(config: StratoMod) -> Tuple[List[RunKeysTrain], List[RunKeysTest]]:
    models = [
        ((model_key, filter_key, data_key), rest)
        for model_key, model in config.ebm_runs.items()
        for filter_key in model.filter
        for data_key, rest in model.inputs.items()
    ]
    train = [
        RunKeysTrain(*meta, train_key)
        for (meta, rest) in models
        for train_key in rest.train
    ]
    test = [
        RunKeysTest(*meta, test_key, test.input_key)
        for (meta, rest) in models
        for test_key, test in rest.test.items()
    ]
    return (train, test)


def test_has_bench(config: StratoMod, runs: RunKeysTest) -> bool:
    return config.inputs[runs.input_key].benchmark is not None


def partition_test_set(
    config: StratoMod,
    test_set: List[RunKeysTest],
) -> Tuple[List, List]:
    unlabeled, labeled = partition(lambda t: test_has_bench(config, t), test_set)
    return list(unlabeled), list(labeled)


def all_refset_keys(
    config: StratoMod,
    ks: Sequence[Union[RunKeysTest, RunKeysTrain]],
) -> List[str]:
    return list(map(lambda x: inputkey_to_refsetkey(config, x.input_key), ks))


def all_input_summary_files(
    config: StratoMod,
    labeled_target: InputFiles,
    unlabeled_target: InputFiles,
) -> InputFiles:
    def labeled_targets(target: InputFiles, key_set: List[RunKeysTrain]):
        return expand(
            target,
            zip,
            model_key=map(lambda x: x.model_key, key_set),
            filter_key=map(lambda x: x.filter_key, key_set),
            input_key=map(lambda x: x.input_key, key_set),
            refset_key=all_refset_keys(config, key_set),
        )

    def test_targets(target: InputFiles, key_set: List[RunKeysTest]) -> InputFiles:
        return expand(
            target,
            zip,
            model_key=map(lambda x: x.model_key, key_set),
            filter_key=map(lambda x: x.filter_key, key_set),
            input_key=map(lambda x: x.test_key, key_set),
            refset_key=all_refset_keys(config, key_set),
        )

    # TODO validate that all inputs in train set have benchmark defined
    train_set, test_set = lookup_run_sets(config)
    unlabeled_test_set, labeled_test_set = partition_test_set(config, test_set)

    return (
        labeled_targets(labeled_target, train_set)
        + test_targets(labeled_target, labeled_test_set)
        + test_targets(unlabeled_target, unlabeled_test_set)
    )


def all_ebm_files(
    config: StratoMod,
    train_target: InputFiles,
    labeled_test_target: InputFiles,
    unlabeled_test_target: InputFiles,
) -> InputFiles:
    def test_targets(path: InputFiles, key_set: List[RunKeysTest]) -> InputFiles:
        return expand(
            path,
            zip,
            model_key=map(lambda x: x.model_key, key_set),
            filter_key=map(lambda x: x.filter_key, key_set),
            data_keys=map(lambda x: x.data_key, key_set),
            input_key=map(lambda x: x.input_key, key_set),
            test_key=map(lambda x: x.test_key, key_set),
            refset_key=all_refset_keys(config, train_set),
        )

    train_set, test_set = lookup_run_sets(config)
    unlabeled_test_set, labeled_test_set = partition_test_set(config, test_set)
    train = expand(
        train_target,
        zip,
        model_key=map(lambda x: x.model_key, train_set),
        filter_key=map(lambda x: x.filter_key, train_set),
        data_key=map(lambda x: x.data_key, train_set),
    )

    # TODO these should eventually point to the test summary htmls
    labeled_test = test_targets(
        labeled_test_target,
        labeled_test_set,
    )
    unlabeled_test = test_targets(
        unlabeled_test_target,
        unlabeled_test_set,
    )

    return train + labeled_test + unlabeled_test


# ------------------------------------------------------------------------------
# global lookup


# def walk_dict(d: StratoMod, keys: List[str]) -> Any:
#     return d if len(keys) == 0 else walk_dict(d[keys[0]], keys[1:])


# def lookup_config(config: StratoMod, keys: List[str]) -> Any:
#     return walk_dict(config, keys)


def chr_index_to_str(i: int) -> str:
    return "X" if i == 23 else ("Y" if i == 24 else str(i))


def chr_indices_to_name(prefix: str, xs: Set[int]) -> Set[str]:
    return set(f"{prefix}{chr_index_to_str(i)}" for i in xs)


# ------------------------------------------------------------------------------
# lookup from ref_key


# def refkey_to_ref(config: StratoMod, args: List[str], ref_key: str) -> Any:
#     return config, ["references", ref_key, *args])


# def refkey_to_benchmark(
#     config: StratoMod,
#     args: List[str],
#     bench_key: str,
#     ref_key: str,
# ) -> Any:
#     return config.references[ref_key].benchmarks[bench_key]


# ------------------------------------------------------------------------------
# lookup from refset_key


# def refsetkey_to_benchmark(
#     config: StratoMod,
#     args: List[str],
#     bench_key: str,
#     refset_key: str,
# ) -> Any:
#     return compose(
#         partial(refkey_to_benchmark, config, args, bench_key),
#         partial(refsetkey_to_refkey, config),
#     )(refset_key)


def refsetkey_to_refset(config: StratoMod, args: List[str], key: RefsetKey) -> RefSet:
    return config.reference_sets[key]


def refsetkey_to_refkey(config: StratoMod, key: RefsetKey) -> RefKey:
    return config.reference_sets[key].ref


def refsetkey_to_ref(config: StratoMod, key: RefsetKey) -> Reference:
    return config.references[refsetkey_to_refkey(config, key)]


def refsetkey_to_chr_indices(config: StratoMod, key: RefsetKey) -> Set[int]:
    f = config.reference_sets[key].chr_filter
    return set(range(1, 25)) if len(f) == 0 else f


# def refsetkey_to_chr_prefix(config: StratoMod, args: List[str], key: RefsetKey) -> str:
#     ref_key = refsetkey_to_refkey(config, key)
#     return refkey_to_ref(config, [*args, "chr_prefix"], ref_key)


def refsetkey_to_sdf_chr_filter(config: StratoMod, key: RefsetKey) -> Set[str]:
    indices = refsetkey_to_chr_indices(config, key)
    prefix = refsetkey_to_ref(config, key).sdf.chr_prefix
    return chr_indices_to_name(prefix, indices)


# ------------------------------------------------------------------------------
# lookup from input_key

# This is tricky because 99% of the time I want to think of train and test
# inputs as being part of a flat list; therefore, "input_key" means "train_key
# or test_key" (most of the time)


# def inputkey_to_input(config: StratoMod, args: List[str], input_key: str) -> VCFInput:
def inputkey_to_input(config: StratoMod, input_key: InputKey) -> VCFInput:
    return config.inputs[input_key]


# def inputkey_to_shared(config: StratoMod, input_key: str, shared_key: str) -> Any:
#     """
#     Given an input key, return a value indicated by 'shared_key', which is a
#     key common to all input dictionaries.
#     """
#     return dict(
#         flatten(
#             [
#                 (train_key, v[shared_key]),
#                 *[(test_key, v[shared_key]) for test_key in v["test"]],
#             ]
#             for train_key, v in config.inputs.items()
#         )
#     )[input_key]


def inputkey_to_refkey(config: StratoMod, input_key: str) -> str:
    return compose(
        partial(refsetkey_to_refkey, config),
        partial(inputkey_to_refsetkey, config),
    )(input_key)


def inputkey_to_refsetkey(config: StratoMod, input_key: InputKey) -> RefsetKey:
    return config.inputs[input_key].refset


# def lookup_inputs(config: StratoMod) -> StratoMod:
#     return config.inputs


# def lookup_train(config: StratoMod, train_key: str) -> StratoMod:
#     return config.inputs.train_key


# def lookup_all_train(config: StratoMod) -> List[Tuple[str, StratoMod]]:
#     return [(k, v["train"]) for k, v in lookup_inputs(config).items()]


# def lookup_all_test(config: StratoMod) -> List[Tuple[str, StratoMod]]:
#     return [ts for i in lookup_inputs(config).values() for ts in i["test"].items()]


# def lookup_test(config: StratoMod, test_key: str) -> StratoMod:
#     return dict(lookup_all_test(config))[test_key]


# def input_train_keys(config: StratoMod):
#     return [*lookup_inputs(config)]


# def input_test_keys(config: StratoMod):
#     return [i[0] for i in lookup_all_test(config)]


def all_refsetkeys(config: StratoMod) -> Set[RefsetKey]:
    return set(v.refset for v in config.inputs.values())


def all_refkeys(config: StratoMod) -> Set[str]:
    return set(map(partial(refsetkey_to_refkey, config), all_refsetkeys(config)))


# TODO this isn't a well defined function in terms of types since the test
# dictionary has several optional fields and the train does not.
# def flat_inputs(config: StratoMod) -> StratoMod:
#     """Return a dictionary of all inputs in the config."""
#     return dict(
#         flatten(
#             [(k, v["train"]), *v["test"].items()] for k, v in config["inputs"].items()
#         )
#     )


def inputkey_to_chr_prefix(config: StratoMod, key: InputKey) -> str:
    return inputkey_to_input(config, key).chr_prefix


def inputkey_to_chr_filter(config: StratoMod, input_key: InputKey) -> Set[int]:
    refset_key = inputkey_to_refsetkey(config, input_key)
    return refsetkey_to_chr_indices(config, refset_key)


def inputkey_to_bench_correction(
    config: StratoMod,
    key: InputKey,
) -> Optional[BenchmarkCorrections]:
    bkey = config.inputs[key].benchmark
    if bkey is not None:
        rkey = inputkey_to_refsetkey(config, key)
        return refsetkey_to_ref(config, rkey).benchmarks[bkey].corrections
    return None


# ------------------------------------------------------------------------------
# variable key lookup


def trainkey_to_variables(config: StratoMod, key: InputKey) -> Dict[VarKey, str]:
    return config.inputs[key].variables


# def testkey_to_variables(config: StratoMod, key: str) -> Dict[str, str]:
#     test = config["inputs"]["test"][key]
#     if "variables" in test:
#         return test["variables"]
#     else:
#         return {}


def runkey_to_variables(
    config: StratoMod,
    mkey: ModelKey,
    rkey: RunKey,
    tkey: TestKey,
) -> Dict[VarKey, str]:
    test = config.ebm_runs[mkey].inputs[rkey].test[tkey]
    default_vars = trainkey_to_variables(config, test.input_key)
    return {**default_vars, **test.variables}


# ------------------------------------------------------------------------------
# ebm run lookup


def lookup_ebm_run(config: StratoMod, run_key: ModelKey) -> Any:
    return config.ebm_runs[run_key]


def ebm_run_train_keys(ebm_run: StratoMod) -> List[str]:
    return [*flatten([[*i] for i in ebm_run.inputs])]


def ebm_run_test_keys(ebm_run: EBMRun) -> List[TestKey]:
    return [
        *flatten([[*ts] for k, v in ebm_run.inputs.items() for ts in v.test.values()])
    ]


# ------------------------------------------------------------------------------
# feature naming lookup (enter if you dare)

# NOTE: The reason this is so convoluted is because I wanted to have a way to
# test if a given feature set is valid, which means we need to know a priori
# what the total possible feature set is. If this seems unjustified, imagine
# running a gigantic configuration on a cluster, only to have it fail at the
# training step after several hours because a feature was named incorrectly.


def lookup_raw_index(config: StratoMod) -> str:
    return config.feature_meta.raw_index


def lookup_bed_cols(config: StratoMod) -> BedIndex:
    return config.feature_meta.bed_index


def bed_cols_ordered(bed_cols: BedIndex) -> List[str]:
    return [bed_cols.chr, bed_cols.start, bed_cols.end]


def bed_cols_indexed(indices: List[int], bed_cols: BedIndex) -> Dict[int, str]:
    return dict(zip(indices, bed_cols_ordered(bed_cols)))


def lookup_bed_cols_ordered(config: StratoMod) -> List[str]:
    return bed_cols_ordered(lookup_bed_cols(config))


def lookup_all_index_cols(config: StratoMod) -> List[str]:
    return [lookup_raw_index(config), *bed_cols_ordered(lookup_bed_cols(config))]


def fmt_feature(prefix: str, rest: str) -> FeatureKey:
    return FeatureKey(f"{prefix}_{rest}")


def fmt_vcf_feature(config: StratoMod, which: str) -> FeatureKey:
    fconf = config.feature_meta.vcf
    return fmt_feature(fconf.prefix, fconf.columns.dict()[which])


def vcf_feature_names(config: StratoMod) -> List[FeatureKey]:
    return [
        *map(
            lambda f: fmt_vcf_feature(config, f),
            config.feature_meta.vcf.columns.dict(),
        )
    ]


def fmt_mappability_feature(config: StratoMod, which: str) -> FeatureKey:
    fconf = config.feature_meta.mappability
    return fmt_feature(fconf.prefix, fconf.suffixes.dict()[which])


def mappability_feature_names(config: StratoMod) -> List[FeatureKey]:
    return [
        *map(
            lambda f: fmt_mappability_feature(config, f),
            config.feature_meta.mappability.suffixes.dict(),
        ),
    ]


def fmt_homopolymer_feature(config: StratoMod, bases: Base, which: str) -> FeatureKey:
    fconf = config.feature_meta.homopolymers
    return fmt_feature(fconf.prefix, f"{bases}_{fconf.suffixes.dict()[which]}")


def homopolymer_feature_names(config: StratoMod) -> List[FeatureKey]:
    fconf = config.feature_meta.homopolymers
    return [
        fmt_homopolymer_feature(config, b, s)
        for (b, s) in product(fconf.bases, list(fconf.suffixes.dict()))
    ]


def fmt_repeat_masker_feature(
    config: StratoMod,
    grp: str,
    fam: Optional[str] = None,
) -> FeatureKey:
    fconf = config.feature_meta.repeat_masker
    rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
    # TODO this is hardcoded for now
    suffix = fconf.suffixes.len
    return fmt_feature(fconf.prefix, f"{rest}_{suffix}")


def repeat_masker_feature_names(config: StratoMod) -> List[FeatureKey]:
    fconf = config.feature_meta.repeat_masker
    fmt = partial(fmt_repeat_masker_feature, config)
    return [
        *[fmt(c) for c in fconf.classes.dict()],
        *[fmt(c, f) for c, fs in fconf.classes.dict().items() for f in fs],
    ]


def fmt_count_feature(prefix: str) -> FeatureKey:
    return fmt_feature(prefix, "count")


def fmt_merged_feature(prefix: str, middle: str, op: BedMergeOp) -> FeatureKey:
    return fmt_feature(prefix, f"{middle}_{op}")


def merged_feature_names(
    prefix: str, names: List[str], ops: Set[BedMergeOp]
) -> List[FeatureKey]:
    return [
        *[fmt_merged_feature(prefix, n, o) for n, o in product(names, ops)],
        fmt_count_feature(prefix),
    ]


def segdup_feature_names(config: StratoMod) -> List[FeatureKey]:
    fconf = config.feature_meta.segdups
    return merged_feature_names(
        fconf.prefix,
        list(fconf.columns.dict().values()),
        fconf.operations,
    )


def fmt_tandem_repeat_base(config: StratoMod, bases: str) -> FeatureKey:
    bs_prefix = config.feature_meta.tandem_repeats.bases_prefix
    return FeatureKey(f"{bs_prefix}_{bases}")


def tandem_repeat_feature_names(config: StratoMod) -> List[FeatureKey]:
    fconf = config.feature_meta.tandem_repeats
    prefix = fconf.prefix
    # TODO weirdly hardcoded in several places
    bases = ["A", "T", "G", "C", "AT", "GC"]
    bs = [fmt_tandem_repeat_base(config, b) for b in bases]
    cs = fconf.columns.dict().values()
    return [
        *merged_feature_names(prefix, [*cs, *bs], fconf.operations),
        *[fmt_feature(prefix, o) for o in fconf.other.dict().values()],
    ]


def all_feature_names(config: StratoMod) -> List[FeatureKey]:
    return [
        *vcf_feature_names(config),
        *mappability_feature_names(config),
        *homopolymer_feature_names(config),
        *repeat_masker_feature_names(config),
        *segdup_feature_names(config),
        *tandem_repeat_feature_names(config),
    ]
