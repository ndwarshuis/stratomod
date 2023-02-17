import re
from typing import (
    TypeVar,
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
    cast,
)
from typing_extensions import Annotated
from more_itertools import flatten, duplicates_everseen, unzip, partition
from itertools import product
from functools import partial
from .functional import maybe, compose
from snakemake.io import expand, InputFiles  # type: ignore
from pydantic import BaseModel as PydanticBaseModel
from pydantic import (
    validator,
    HttpUrl,
    PositiveFloat,
    NonNegativeInt,
    conint,
    constr,
    conset,
    confloat,
)
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


class BaseModel(PydanticBaseModel):
    class Config:
        validate_all = True
        extra = "forbid"
        frozen = True


class Paths(BaseModel):
    resources: str
    results: str


ChrIndex = Annotated[int, conint(ge=0, le=24)]


class RefSet(BaseModel):
    ref: RefKey
    chr_filter: Set[ChrIndex]


class CatVar(BaseModel):
    levels: Annotated[Set[str], conset(str, min_items=1)]


class Truncation(BaseModel):
    lower: Optional[float] = None
    upper: Optional[float] = None

    @validator("upper")
    def lower_less_than_upper(cls, v: Optional[float], values) -> Optional[float]:
        lower = values["lower"]
        if v is not None and values["lower"] is not None:
            assert lower <= v
        return v

    def in_range(self, x: float) -> bool:
        lower = self.lower
        upper = self.upper
        return (True if lower is None else lower <= x) and (
            True if upper is None else x <= upper
        )


class ContVar(Truncation):
    lower: Optional[float] = None
    upper: Optional[float] = None


class Variables(BaseModel):
    continuous: Dict[VarKey, ContVar]
    categorical: Dict[VarKey, CatVar]

    def validate_variable(self, varname: VarKey, varval: str) -> None:
        cats = self.categorical
        conts = self.continuous
        if varname in cats:
            assert (
                varval in cats[varname].levels
            ), f"'{varname}' has an invalid level '{varval}'"
        elif varname in conts:
            assert varval.isnumeric(), f"'{varname}' not a number"
            assert conts[varname].in_range(float(varval))
        else:
            assert False, f"'{varname}' not a valid variable name"


class FormatFields(BaseModel):
    vaf: Optional[str]
    dp: Optional[str]
    gt: Optional[str]
    gq: Optional[str]


class VCFInput(BaseModel):
    refset: RefsetKey
    chr_prefix: str
    url: Optional[HttpUrl]
    benchmark: Optional[BenchKey]
    variables: Dict[VarKey, str]
    format_fields: Optional[FormatFields]
    max_ref: Annotated[int, conint(ge=0)] = 50
    max_alt: Annotated[int, conint(ge=0)] = 50


class TestDataInput(BaseModel):
    input_key: InputKey
    variables: Dict[VarKey, str]


class ModelRun(BaseModel):
    train: Annotated[Set[InputKey], conset(InputKey, min_items=1)]
    test: Dict[TestKey, TestDataInput]


Fraction = Annotated[float, confloat(ge=0, le=1, allow_inf_nan=False)]


class EBMMiscParams(BaseModel):
    downsample: Optional[Fraction] = None


class EBMSplitParams(BaseModel):
    test_size: Fraction = 0.2
    random_state: Optional[Fraction] = None


# TODO these all need constraints
class EBMClassifierParams(BaseModel):
    max_bins: NonNegativeInt = 256
    max_interaction_bins: NonNegativeInt = 32
    binning: Binning = Binning.QUANTILE
    outer_bags: NonNegativeInt = 8
    inner_bags: NonNegativeInt = 0
    learning_rate: PositiveFloat = 0.01
    validation_size: PositiveFloat = 0.15
    early_stopping_rounds: NonNegativeInt = 50
    early_stopping_tolerance: PositiveFloat = 0.0001
    max_rounds: NonNegativeInt = 5000
    min_samples_leaf: NonNegativeInt = 2
    max_leaves: NonNegativeInt = 3
    random_state: Optional[int] = None


class EBMSettings(BaseModel):
    misc_parameters: EBMMiscParams = EBMMiscParams()
    split_parameters: EBMSplitParams = EBMSplitParams()
    classifier_parameters: EBMClassifierParams = EBMClassifierParams()


class Visualization(BaseModel):
    truncate: Truncation = Truncation()
    plot_type: PlotType = PlotType.STEP
    split_missing: Optional[Fraction] = None


class Feature(BaseModel):
    feature_type: FeatureType = FeatureType.CONTINUOUS
    fill_na: Optional[float] = 0.0
    alt_name: Optional[FeatureKey] = None
    visualization: Visualization = Visualization()
    transform: Optional[Transform] = None


class FeaturePair(BaseModel):
    f1: FeatureKey
    f2: FeatureKey


class Model(BaseModel):
    runs: Dict[RunKey, ModelRun]
    filter: Set[FilterKey]
    ebm_settings: EBMSettings
    # TODO only allow actual error labels here
    error_labels: Annotated[Set[Label], conset(Label, min_items=1)]
    filtered_are_candidates: bool
    interactions: Union[NonNegativeInt, Set[Union[FeatureKey, FeaturePair]]] = 0
    features: Dict[FeatureKey, Feature]

    @validator("features")
    def model_has_matching_alt_features(
        cls, fs: Dict[FeatureKey, Feature]
    ) -> Dict[FeatureKey, Feature]:
        for k, v in fs.items():
            alt = v.alt_name
            if alt is not None:
                prefix = assert_match("^[^_]+", k)
                alt_prefix = assert_match("^[^_]+", alt)
                assert alt_prefix == prefix, f"Alt prefix must match for {k}"
        return fs


class BedFile(BaseModel):
    url: Optional[HttpUrl]
    chr_prefix: str


class Strats(BaseModel):
    mhc: BedFile


class BenchmarkCorrections(BaseModel):
    strip_IPS: bool


class Benchmark(BaseModel):
    vcf_url: Optional[HttpUrl]
    bed_url: Optional[HttpUrl]
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
    benchmarks: Dict[BenchKey, Benchmark]


Prefix = Annotated[str, constr(regex="[A-Z]+")]

NonEmptyStr = Annotated[str, constr(min_length=1)]


class BedIndex(BaseModel):
    chr: NonEmptyStr
    start: NonEmptyStr
    end: NonEmptyStr


class VCFColumns(BaseModel):
    input: NonEmptyStr
    qual: NonEmptyStr
    filter: NonEmptyStr
    info: NonEmptyStr
    gt: NonEmptyStr
    gq: NonEmptyStr
    dp: NonEmptyStr
    vaf: NonEmptyStr
    len: NonEmptyStr


class VCFMeta(BaseModel):
    prefix: Prefix
    columns: VCFColumns


class MapSuffixes(BaseModel):
    low: NonEmptyStr
    high: NonEmptyStr


class MapMeta(BaseModel):
    prefix: Prefix
    suffixes: MapSuffixes


class HomopolySuffixes(BaseModel):
    len: NonEmptyStr
    imp_frac: NonEmptyStr


class HomopolyMeta(BaseModel):
    prefix: Prefix
    bases: Set[Base]
    suffixes: HomopolySuffixes


class RMSKSuffixes(BaseModel):
    len: str


class RMSKClasses(BaseModel):
    SINE: Set[NonEmptyStr] = set()
    LINE: Set[NonEmptyStr] = set()
    LTR: Set[NonEmptyStr] = set()
    Satellite: Set[NonEmptyStr] = set()


class RMSKMeta(BaseModel):
    prefix: Prefix
    suffixes: RMSKSuffixes
    classes: RMSKClasses


class SegDupsColumns(BaseModel):
    alignL: NonEmptyStr
    fracMatchIndel: NonEmptyStr


class SegDupsMeta(BaseModel):
    prefix: Prefix
    columns: SegDupsColumns
    operations: Set[BedMergeOp]


class TandemRepeatColumns(BaseModel):
    period: NonEmptyStr
    copyNum: NonEmptyStr
    perMatch: NonEmptyStr
    perIndel: NonEmptyStr
    score: NonEmptyStr


class TandemRepeatOther(BaseModel):
    len: NonEmptyStr


class TandemRepeatMeta(BaseModel):
    prefix: Prefix
    bases_prefix: NonEmptyStr
    operations: Set[BedMergeOp]
    columns: TandemRepeatColumns
    other: TandemRepeatOther


class FeatureMeta(BaseModel):
    label: NonEmptyStr
    raw_index: NonEmptyStr
    bed_index: BedIndex
    vcf: VCFMeta
    mappability: MapMeta
    homopolymers: HomopolyMeta
    repeat_masker: RMSKMeta
    segdups: SegDupsMeta
    tandem_repeats: TandemRepeatMeta


class Tools(BaseModel):
    repseq: HttpUrl


def flatten_features(fs: Dict[FeatureKey, Feature]) -> List[FeatureKey]:
    return [k if v.alt_name is None else v.alt_name for k, v in fs.items()]


class StratoMod(BaseModel):
    paths: Paths
    tools: Tools
    feature_meta: FeatureMeta
    variables: Variables
    references: Dict[RefKey, Reference]
    reference_sets: Dict[RefsetKey, RefSet]
    inputs: Dict[InputKey, VCFInput]
    models: Dict[ModelKey, Model]

    @validator("reference_sets", each_item=True)
    def refsets_have_valid_refkeys(cls, v: RefSet, values) -> RefSet:
        assert (
            v.ref in values["references"]
        ), f"'{v.ref}' does not refer to a valid reference"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_refsetkeys(cls, v: VCFInput, values) -> VCFInput:
        assert (
            v.refset in values["reference_sets"]
        ), f"'{v.refset}' does not refer to a valid reference set"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_benchkeys(cls, v: VCFInput, values) -> VCFInput:
        if v.benchmark is not None:
            ref_key = values["reference_sets"][v.refset].ref
            ref_benchmarks = values["references"][ref_key].benchmarks
            assert (
                v.benchmark in ref_benchmarks
            ), f"'{v.benchmark}' does not refer to a valid benchmark"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_variables(cls, v: VCFInput, values) -> VCFInput:
        var_root = cast(Variables, values["variables"])
        for varname, varval in v.variables.items():
            var_root.validate_variable(varname, varval)
        return v

    @validator("models", each_item=True)
    def models_have_valid_features(cls, v: Model, values) -> Model:
        features = all_feature_names(values["feature_meta"])
        # TODO add variable names here
        assert_subset(set(v.features), features)
        # assert set(v.features) <= features
        return v

    @validator("models", each_item=True)
    def models_have_valid_features_alt(cls, v: Model) -> Model:
        # TODO dry?
        # TODO need to also compare variable names
        assert_no_dups(flatten_features(v.features), "Duplicated features")
        return v

    @validator("models", each_item=True)
    def models_have_valid_interactions(cls, v: Model, values) -> Model:
        if isinstance(v.interactions, set):
            # TODO add variable names here
            features = set(flatten_features(v.features))
            interactions = set(
                flatten(
                    [i.f1, i.f2] if isinstance(i, FeaturePair) else [i]
                    for i in v.interactions
                )
            )
            assert_subset(interactions, features)
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_train(cls, v: Model, values) -> Model:
        train = [t for r in v.runs.values() for t in r.train]
        assert_subset(set(train), set(values["inputs"]))
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_test(cls, v: Model, values) -> Model:
        tests = [t.input_key for r in v.runs.values() for t in r.test.values()]
        assert_subset(set(tests), set(values["inputs"]))
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_test_variables(cls, v: Model, values) -> Model:
        var_root = cast(Variables, values["variables"])
        varpairs = [
            (varname, varval)
            for r in v.runs.values()
            for t in r.test.values()
            for varname, varval in t.variables.items()
        ]
        for varname, varval in varpairs:
            var_root.validate_variable(varname, varval)
        return v


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


X = TypeVar("X")


def assert_subset(xs: Set[X], ys: Set[X]) -> None:
    assert xs <= ys, f"not a subset - extra members: {xs - ys}"


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


def validate_input_keys(
    config: StratoMod,
    xs: List[InputKey],
    validate_bench: bool,
) -> None:
    assert set(xs) <= set(config.inputs)
    validate_refset_keys(config, [inputkey_to_refsetkey(config, x) for x in xs])
    # TODO validate variables
    if validate_bench:
        bench_pairs = [(x, config.inputs[x].benchmark) for x in xs]
        no_bench = [i for i, b in bench_pairs if b is None]
        assert len(no_bench) == 0
        validate_bench_keys(
            config,
            [
                (inputkey_to_refkey(config, input_key), bench_key)
                for input_key, bench_key in bench_pairs
                if bench_key is not None
            ],
        )


def validate_ref_keys(config: StratoMod, xs: List[RefKey]) -> None:
    assert set(xs) <= set(config.references)


def validate_bench_keys(config: StratoMod, xs: List[Tuple[RefKey, BenchKey]]) -> None:
    assert set(xs) <= set(config.references)


def validate_refset_keys(config: StratoMod, xs: List[RefsetKey]) -> None:
    assert set(xs) <= set(config.reference_sets)
    validate_ref_keys(config, [refsetkey_to_refkey(config, x) for x in xs])


# def validate_ebm_test_inputs(xs: List[RunKeysTest], wants_bench: bool) -> None:
#     pass


def lookup_run_sets(config: StratoMod) -> Tuple[List[RunKeysTrain], List[RunKeysTest]]:
    models = [
        ((model_key, filter_key, data_key), rest)
        for model_key, model in config.models.items()
        for filter_key in model.filter
        for data_key, rest in model.runs.items()
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

    # validate_ebm_train_inputs(config, train_set)
    # validate_ebm_test_inputs(labeled_test_set, True)
    # validate_ebm_test_inputs(unlabeled_test_set, False)

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


def inputkey_to_refkey(config: StratoMod, input_key: InputKey) -> RefKey:
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


# def trainkey_to_variables(config: StratoMod, key: InputKey) -> Dict[VarKey, str]:
#     return config.inputs[key].variables


# def testkey_to_variables(config: StratoMod, key: str) -> Dict[str, str]:
#     test = config["inputs"]["test"][key]
#     if "variables" in test:
#         return test["variables"]
#     else:
#         return {}


# def runkey_to_variables(
#     config: StratoMod,
#     mkey: ModelKey,
#     rkey: RunKey,
#     tkey: TestKey,
# ) -> Dict[VarKey, str]:
#     test = config.ebm_runs[mkey].inputs[rkey].test[tkey]
#     default_vars = trainkey_to_variables(config, test.input_key)
#     return {**default_vars, **test.variables}


# ------------------------------------------------------------------------------
# ebm run lookup


def lookup_ebm_run(config: StratoMod, run_key: ModelKey) -> Any:
    return config.models[run_key]


def ebm_run_train_keys(ebm_run: StratoMod) -> List[str]:
    return [*flatten([[*i] for i in ebm_run.inputs])]


def ebm_run_test_keys(ebm_run: Model) -> List[TestKey]:
    return [
        *flatten([[*ts] for k, v in ebm_run.runs.items() for ts in v.test.values()])
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


def fmt_vcf_feature(config: FeatureMeta, which: str) -> FeatureKey:
    fconf = config.vcf
    return fmt_feature(fconf.prefix, fconf.columns.dict()[which])


def vcf_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    return [
        *map(
            lambda f: fmt_vcf_feature(config, f),
            config.vcf.columns.dict(),
        )
    ]


def fmt_mappability_feature(config: FeatureMeta, which: str) -> FeatureKey:
    fconf = config.mappability
    return fmt_feature(fconf.prefix, fconf.suffixes.dict()[which])


def mappability_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    return [
        *map(
            lambda f: fmt_mappability_feature(config, f),
            config.mappability.suffixes.dict(),
        ),
    ]


def fmt_homopolymer_feature(config: FeatureMeta, bases: Base, which: str) -> FeatureKey:
    fconf = config.homopolymers
    return fmt_feature(fconf.prefix, f"{bases.value}_{fconf.suffixes.dict()[which]}")


def homopolymer_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    fconf = config.homopolymers
    return [
        fmt_homopolymer_feature(config, b, s)
        for (b, s) in product(fconf.bases, list(fconf.suffixes.dict()))
    ]


def fmt_repeat_masker_feature(
    fconf: RMSKMeta,
    grp: str,
    fam: Optional[str],
) -> FeatureKey:
    rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
    # TODO this is hardcoded for now
    suffix = fconf.suffixes.len
    return fmt_feature(fconf.prefix, f"{rest}_{suffix}")


def repeat_masker_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    fconf = config.repeat_masker

    def fmt(grp: str, fam: Optional[str]) -> FeatureKey:
        return fmt_repeat_masker_feature(fconf, grp, fam)

    return [
        *[fmt(c, None) for c in fconf.classes.dict()],
        *[fmt(c, f) for c, fs in fconf.classes.dict().items() for f in fs],
    ]


def fmt_count_feature(prefix: str) -> FeatureKey:
    return fmt_feature(prefix, "count")


def fmt_merged_feature(prefix: str, middle: str, op: BedMergeOp) -> FeatureKey:
    return fmt_feature(prefix, f"{middle}_{op.value}")


def merged_feature_names(
    prefix: str, names: List[str], ops: Set[BedMergeOp]
) -> List[FeatureKey]:
    return [
        *[fmt_merged_feature(prefix, n, o) for n, o in product(names, ops)],
        fmt_count_feature(prefix),
    ]


def segdup_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    fconf = config.segdups
    return merged_feature_names(
        fconf.prefix,
        list(fconf.columns.dict().values()),
        fconf.operations,
    )


def fmt_tandem_repeat_base(config: FeatureMeta, bases: str) -> FeatureKey:
    bs_prefix = config.tandem_repeats.bases_prefix
    return FeatureKey(f"{bs_prefix}_{bases}")


def tandem_repeat_feature_names(config: FeatureMeta) -> List[FeatureKey]:
    fconf = config.tandem_repeats
    prefix = fconf.prefix
    # TODO weirdly hardcoded in several places
    bases = ["A", "T", "G", "C", "AT", "GC"]
    bs = [fmt_tandem_repeat_base(config, b) for b in bases]
    cs = fconf.columns.dict().values()
    return [
        *merged_feature_names(prefix, [*cs, *bs], fconf.operations),
        *[fmt_feature(prefix, o) for o in fconf.other.dict().values()],
    ]


def all_feature_names(config: FeatureMeta) -> Set[FeatureKey]:
    return set(
        [
            *vcf_feature_names(config),
            *mappability_feature_names(config),
            *homopolymer_feature_names(config),
            *repeat_masker_feature_names(config),
            *segdup_feature_names(config),
            *tandem_repeat_feature_names(config),
        ]
    )
