import re
from pathlib import Path
from typing import (
    TypeVar,
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
from .functional import maybe
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
    resources: Path
    results: Path


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
        return maybe(True, lambda y: y <= x, self.lower) and maybe(
            True, lambda y: x <= y, self.upper
        )


class ContVar(Truncation):
    lower: Optional[float] = None
    upper: Optional[float] = None


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

    def bed_cols_ordered(self) -> List[str]:
        return [self.chr, self.start, self.end]

    def bed_cols_indexed(self, indices: Tuple[int, int, int]) -> Dict[int, str]:
        return dict(zip(indices, self.bed_cols_ordered()))


class FeatureGroup(BaseModel):
    prefix: Prefix

    def fmt_feature(self, rest: str) -> FeatureKey:
        return FeatureKey(f"{self.prefix}_{rest}")


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


class VCFMeta(FeatureGroup):
    columns: VCFColumns

    @validator("columns")
    def prefix_columns(cls, columns: VCFColumns, values) -> VCFColumns:
        return VCFColumns(
            **{k: f"{values['prefix']}_{v}" for k, v in columns.dict().items()}
        )

    def feature_names(self) -> List[FeatureKey]:
        return list(self.columns.dict().values())


class MapSuffixes(BaseModel):
    low: NonEmptyStr
    high: NonEmptyStr


class MapMeta(FeatureGroup):
    suffixes: MapSuffixes

    @property
    def low(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.low)

    @property
    def high(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.high)

    def feature_names(self) -> List[FeatureKey]:
        return [self.low, self.high]


class HomopolySuffixes(BaseModel):
    len: NonEmptyStr
    imp_frac: NonEmptyStr


class HomopolyMeta(FeatureGroup):
    bases: Set[Base]
    suffixes: HomopolySuffixes

    # TODO weakly typed
    def fmt_name(self, bases: Base, which: str) -> FeatureKey:
        return self.fmt_feature(f"{bases.value}_{self.suffixes.dict()[which]}")

    def feature_names(self) -> List[FeatureKey]:
        return [
            self.fmt_name(b, s)
            for (b, s) in product(self.bases, list(self.suffixes.dict()))
        ]


class RMSKSuffixes(BaseModel):
    len: str


class RMSKClasses(BaseModel):
    SINE: Set[NonEmptyStr] = set()
    LINE: Set[NonEmptyStr] = set()
    LTR: Set[NonEmptyStr] = set()
    Satellite: Set[NonEmptyStr] = set()


class RMSKMeta(FeatureGroup):
    suffixes: RMSKSuffixes
    classes: RMSKClasses

    # TODO weakly typed
    def fmt_name(
        self,
        grp: str,
        fam: Optional[str],
    ) -> FeatureKey:
        rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
        # TODO this is hardcoded for now
        suffix = self.suffixes.len
        return self.fmt_feature(f"{rest}_{suffix}")

    def feature_names(self) -> List[FeatureKey]:
        def fmt(grp: str, fam: Optional[str]) -> FeatureKey:
            return self.fmt_name(grp, fam)

        return [
            *[fmt(c, None) for c in self.classes.dict()],
            *[fmt(c, f) for c, fs in self.classes.dict().items() for f in fs],
        ]


class SegDupsColumns(BaseModel):
    alignL: NonEmptyStr
    fracMatchIndel: NonEmptyStr


class MergedFeatureGroup(FeatureGroup):
    operations: Set[BedMergeOp]

    def fmt_count_feature(self) -> FeatureKey:
        return self.fmt_feature("count")

    def fmt_merged_feature(self, middle: str, op: BedMergeOp) -> FeatureKey:
        return self.fmt_feature(f"{middle}_{op.value}")

    def merged_feature_names(self, names: List[str]) -> List[FeatureKey]:
        return [
            *[
                self.fmt_merged_feature(n, o)
                for n, o in product(names, self.operations)
            ],
            self.fmt_count_feature(),
        ]


class SegDupsMeta(MergedFeatureGroup):
    columns: SegDupsColumns

    def feature_names(self) -> List[FeatureKey]:
        return self.merged_feature_names(list(self.columns.dict().values()))


class TandemRepeatColumns(BaseModel):
    period: NonEmptyStr
    copyNum: NonEmptyStr
    perMatch: NonEmptyStr
    perIndel: NonEmptyStr
    score: NonEmptyStr


class TandemRepeatOther(BaseModel):
    len: NonEmptyStr


class TandemRepeatMeta(MergedFeatureGroup):
    bases_prefix: NonEmptyStr
    columns: TandemRepeatColumns
    other: TandemRepeatOther

    def fmt_name_base(self, bases: str) -> FeatureKey:
        bs_prefix = self.bases_prefix
        return FeatureKey(f"{bs_prefix}_{bases}")

    def feature_names(self) -> List[FeatureKey]:
        # TODO weirdly hardcoded in several places
        bases = ["A", "T", "G", "C", "AT", "GC"]
        bs = [self.fmt_name_base(b) for b in bases]
        cs = self.columns.dict().values()
        return [
            *self.merged_feature_names([*cs, *bs]),
            *[self.fmt_feature(o) for o in self.other.dict().values()],
        ]


class Variables(FeatureGroup):
    continuous: Dict[VarKey, ContVar]
    categorical: Dict[VarKey, CatVar]

    def all_keys(self) -> List[VarKey]:
        return list(self.continuous) + list(self.categorical)

    def feature_names(self) -> List[FeatureKey]:
        return [self.fmt_feature(x) for x in self.all_keys()]

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
    variables: Variables

    def all_index_cols(self) -> List[str]:
        return [self.raw_index, *self.bed_index.bed_cols_ordered()]

    def all_feature_names(self) -> Set[FeatureKey]:
        return set(
            [
                *self.vcf.feature_names(),
                *self.mappability.feature_names(),
                *self.homopolymers.feature_names(),
                *self.repeat_masker.feature_names(),
                *self.segdups.feature_names(),
                *self.tandem_repeats.feature_names(),
                *self.variables.feature_names(),
            ]
        )


class Tools(BaseModel):
    repseq: HttpUrl


def flatten_features(fs: Dict[FeatureKey, Feature]) -> List[FeatureKey]:
    return [k if v.alt_name is None else v.alt_name for k, v in fs.items()]


class StratoMod(BaseModel):
    paths: Paths
    tools: Tools
    feature_meta: FeatureMeta
    references: Dict[RefKey, Reference]
    reference_sets: Dict[RefsetKey, RefSet]
    inputs: Dict[InputKey, VCFInput]
    models: Dict[ModelKey, Model]

    @validator("reference_sets", each_item=True)
    def refsets_have_valid_refkeys(cls, v: RefSet, values) -> RefSet:
        if "reference" in values:
            assert (
                v.ref in values["references"]
            ), f"'{v.ref}' does not refer to a valid reference"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_refsetkeys(cls, v: VCFInput, values) -> VCFInput:
        if "reference_sets" in values:
            assert (
                v.refset in values["reference_sets"]
            ), f"'{v.refset}' does not refer to a valid reference set"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_benchkeys(cls, v: VCFInput, values) -> VCFInput:
        if (
            v.benchmark is not None
            and "reference_sets" in values
            and "references" in values
        ):
            ref_key = values["reference_sets"][v.refset].ref
            ref_benchmarks = values["references"][ref_key].benchmarks
            assert (
                v.benchmark in ref_benchmarks
            ), f"'{v.benchmark}' does not refer to a valid benchmark"
        return v

    @validator("inputs", each_item=True)
    def inputs_have_valid_variables(cls, v: VCFInput, values) -> VCFInput:
        if "feature_meta" in values:
            var_root = cast(FeatureMeta, values["feature_meta"]).variables
            for varname, varval in v.variables.items():
                var_root.validate_variable(varname, varval)
        return v

    @validator("models", each_item=True)
    def models_have_valid_features(cls, v: Model, values) -> Model:
        if "feature_meta" in values:
            features = values["feature_meta"].all_feature_names()
            assert_subset(set(v.features), features)
        return v

    @validator("models", each_item=True)
    def models_have_valid_features_alt(cls, v: Model) -> Model:
        # TODO dry?
        assert_no_dups(flatten_features(v.features), "Duplicated features")
        return v

    @validator("models", each_item=True)
    def models_have_valid_interactions(cls, v: Model, values) -> Model:
        if isinstance(v.interactions, set):
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
        if "inputs" in values:
            train = [t for r in v.runs.values() for t in r.train]
            assert_subset(set(train), set(values["inputs"]))
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_test(cls, v: Model, values) -> Model:
        if "inputs" in values:
            tests = [t.input_key for r in v.runs.values() for t in r.test.values()]
            assert_subset(set(tests), set(values["inputs"]))
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_test_variables(cls, v: Model, values) -> Model:
        if "feature_meta" in values:
            var_root = cast(FeatureMeta, values["feature_meta"]).variables
            varpairs = [
                (varname, varval)
                for r in v.runs.values()
                for t in r.test.values()
                for varname, varval in t.variables.items()
            ]
            assert_subset(set([p[0] for p in varpairs]), set(var_root.all_keys()))
            for varname, varval in varpairs:
                var_root.validate_variable(varname, varval)
        return v

    # TODO might be nice to alert user when they try and test a variable that's
    # not included in the train set, since this will likely screw with things
    # in weird ways

    def refsetkey_to_ref(self, key: RefsetKey) -> Reference:
        return self.references[self.refsetkey_to_refkey(key)]

    def refsetkey_to_refset(self, key: RefsetKey) -> RefSet:
        return self.reference_sets[key]

    def inputkey_to_input(self, input_key: InputKey) -> VCFInput:
        return self.inputs[input_key]

    def refsetkey_to_refkey(self, key: RefsetKey) -> RefKey:
        return self.refsetkey_to_refset(key).ref

    def refsetkey_to_chr_indices(self, key: RefsetKey) -> Set[int]:
        f = self.refsetkey_to_refset(key).chr_filter
        return set(range(1, 25)) if len(f) == 0 else f

    def refsetkey_to_sdf_chr_filter(self, key: RefsetKey) -> Set[str]:
        indices = self.refsetkey_to_chr_indices(key)
        prefix = self.refsetkey_to_ref(key).sdf.chr_prefix
        return chr_indices_to_name(prefix, indices)

    def inputkey_to_refkey(self, key: InputKey) -> RefKey:
        rk = self.inputkey_to_refsetkey(key)
        return self.refsetkey_to_refkey(rk)

    def inputkey_to_refsetkey(self, key: InputKey) -> RefsetKey:
        return self.inputkey_to_input(key).refset

    def all_refsetkeys(self) -> Set[RefsetKey]:
        return set(v.refset for v in self.inputs.values())

    def all_refkeys(self) -> Set[str]:
        return set(map(self.refsetkey_to_refkey, self.all_refsetkeys()))

    def inputkey_to_chr_prefix(self, key: InputKey) -> str:
        return self.inputkey_to_input(key).chr_prefix

    def inputkey_to_chr_filter(self, input_key: InputKey) -> Set[int]:
        refset_key = self.inputkey_to_refsetkey(input_key)
        return self.refsetkey_to_chr_indices(refset_key)

    def inputkey_to_variables(self, input_key: InputKey) -> Dict[VarKey, str]:
        return self.inputkey_to_input(input_key).variables

    def inputkey_to_bench_correction(
        self,
        key: InputKey,
    ) -> Optional[BenchmarkCorrections]:
        bkey = self.inputkey_to_input(key).benchmark
        if bkey is not None:
            rkey = self.inputkey_to_refsetkey(key)
            return self.refsetkey_to_ref(rkey).benchmarks[bkey].corrections
        return None

    def testkey_to_variables(
        self,
        mkey: ModelKey,
        rkey: RunKey,
        tkey: TestKey,
    ) -> Dict[VarKey, str]:
        test = self.models[mkey].runs[rkey].test[tkey]
        return {**test.variables, **self.inputkey_to_variables(test.input_key)}


# ------------------------------------------------------------------------------
# too useful...


def fmt_strs(ss: Iterable[str]) -> str:
    return ", ".join(ss)


# ------------------------------------------------------------------------------
# assertions


X = TypeVar("X")


def assert_empty(xs: Collection[X], msg: str) -> None:
    assert len(xs) == 0, f"{msg}: {xs}"


def assert_no_dups(xs: Collection[X], msg: str) -> None:
    assert_empty(set(duplicates_everseen(xs)), msg)


def assert_match(pat: str, s: str) -> str:
    res = re.match(pat, s)
    assert res is not None, f"match failed for pattern '{pat}' and query '{s}'"
    return res[0]


def assert_subset(xs: Set[X], ys: Set[X]) -> None:
    assert xs <= ys, f"not a subset - extra members: {xs - ys}"


# ------------------------------------------------------------------------------
# resources


def attempt_mem_gb(mem_gb: int) -> Callable[[Dict[str, str], int], int]:
    # double initial memory on each attempt
    return lambda wildcards, attempt: mem_gb * 1000 * 2 ** (attempt - 1)


# ------------------------------------------------------------------------------
# expanding targets


def all_benchkeys(config: StratoMod, target: InputFiles) -> InputFiles:
    rs, bs = unzip(
        (config.inputkey_to_refkey(k), b)
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
) -> Tuple[List[RunKeysTest], List[RunKeysTest]]:
    unlabeled, labeled = partition(lambda t: test_has_bench(config, t), test_set)
    return list(unlabeled), list(labeled)


def all_refset_keys(
    config: StratoMod,
    ks: Sequence[Union[RunKeysTest, RunKeysTrain]],
) -> List[str]:
    return list(map(lambda x: config.inputkey_to_refsetkey(x.input_key), ks))


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


def chr_index_to_str(i: ChrIndex) -> str:
    return "X" if i == 23 else ("Y" if i == 24 else str(i))


def chr_indices_to_name(prefix: str, xs: Set[int]) -> Set[str]:
    return set(f"{prefix}{chr_index_to_str(i)}" for i in xs)


# def lookup_ebm_run(config: StratoMod, run_key: ModelKey) -> Any:
#     return config.models[run_key]


# def ebm_run_train_keys(ebm_run: StratoMod) -> List[str]:
#     return [*flatten([[*i] for i in ebm_run.inputs])]


# def ebm_run_test_keys(ebm_run: Model) -> List[TestKey]:
#     return [
#         *flatten([[*ts] for k, v in ebm_run.runs.items() for ts in v.test.values()])
#     ]
