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
# Querykey = NewType("Querykey", str)
UnlabeledQueryKey = NewType("UnlabeledQueryKey", str)
LabeledQueryKey = NewType("LabeledQueryKey", str)
BenchKey = NewType("BenchKey", str)
ModelKey = NewType("ModelKey", str)
RunKey = NewType("RunKey", str)
FeatureKey = NewType("FeatureKey", str)
VarKey = NewType("VarKey", str)

QueryKey = Union[UnlabeledQueryKey, LabeledQueryKey]


class ListEnum(Enum):
    @classmethod
    def all(cls):
        return [x.value for x in cls]


class Base(ListEnum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"


class FilterKey(ListEnum):
    SNV = "SNV"
    INDEL = "INDEL"


class Label(ListEnum):
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


def alternate_constraint(xs: List[str]) -> str:
    return f"({'|'.join(xs)})"


_constraints = {
    # corresponds to a genome reference
    "ref_key": "[^/]+",
    # corresponds to a reference set (reference + chromosome filter + etc)
    "refset_key": "[^/]+",
    # corresponds to a testing vcf
    "test_key": "[^/]+",
    "query_key": "[^/]+",
    "ul_query_key": "[^/]+",
    "l_query_key": "[^/]+",
    # refers to a collection of input data input to an ebm model configuration
    # (composed of multiple train/test keys + associated data)
    "run_key": "[^/]+",
    # refers to a benchmark vcf (within the context of a given reference)
    "bench_key": "[^/]+",
    # refers to an EBM model and its parameters
    "model_key": "[^/]+",
    # refers to the variant type (SNP or INDEL, for now)
    # TODO ...why "filter"? (I can't think of anything better)
    "filter_key": alternate_constraint(FilterKey.all()),
    # refers to a variant benchmarking label (tp, fp, etc)
    "label": alternate_constraint(Label.all()),
    # refers to a nucleotide base
    "base": alternate_constraint(Base.all()),
}

all_wildcards = {k: f"{{{k},{v}}}" for k, v in _constraints.items()}


def wildcard_ext(key, ext):
    return f"{all_wildcards[key]}.{ext}"


def wildcard_format(format_str, *keys):
    return format_str.format(*[all_wildcards[k] for k in keys])


def wildcard_format_ext(format_str, keys, ext):
    return wildcard_format(f"{format_str}.{ext}", *keys)


class BaseModel(PydanticBaseModel):
    class Config:
        validate_all = True
        extra = "forbid"
        frozen = True


class Paths(BaseModel):
    resources: Path
    results: Path
    log: Path


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


class UnlabeledVCFInput(BaseModel):
    refset: RefsetKey
    chr_prefix: str
    url: Optional[HttpUrl]
    variables: Dict[VarKey, str]
    format_fields: Optional[FormatFields]
    max_ref: Annotated[int, conint(ge=0)] = 50
    max_alt: Annotated[int, conint(ge=0)] = 50


class LabeledVCFInput(UnlabeledVCFInput):
    benchmark: BenchKey


VCFInput = Union[UnlabeledVCFInput, LabeledVCFInput]


class TestDataInput(BaseModel):
    query_key: QueryKey
    variables: Dict[VarKey, str]


class ModelRun(BaseModel):
    train: Annotated[Set[LabeledQueryKey], conset(LabeledQueryKey, min_items=1)]
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

    # TODO weakly typed
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


LabeledQueries = Dict[LabeledQueryKey, LabeledVCFInput]
UnlabeledQueries = Dict[UnlabeledQueryKey, UnlabeledVCFInput]


class StratoMod(BaseModel):
    paths: Paths
    tools: Tools
    feature_meta: FeatureMeta
    references: Dict[RefKey, Reference]
    reference_sets: Dict[RefsetKey, RefSet]
    labeled_queries: LabeledQueries
    unlabeled_queries: UnlabeledQueries
    models: Dict[ModelKey, Model]

    @validator("reference_sets", each_item=True)
    def refsets_have_valid_refkeys(cls, v: RefSet, values) -> RefSet:
        if "reference" in values:
            assert (
                v.ref in values["references"]
            ), f"'{v.ref}' does not refer to a valid reference"
        return v

    @validator("labeled_queries", "unlabeled_queries", each_item=True)
    def inputs_have_valid_refsetkeys(cls, v: VCFInput, values) -> VCFInput:
        if "reference_sets" in values:
            assert (
                v.refset in values["reference_sets"]
            ), f"'{v.refset}' does not refer to a valid reference set"
        return v

    @validator("labeled_queries", "unlabeled_queries", each_item=True)
    def inputs_have_valid_variables(cls, v: VCFInput, values) -> VCFInput:
        if "feature_meta" in values:
            var_root = cast(FeatureMeta, values["feature_meta"]).variables
            for varname, varval in v.variables.items():
                var_root.validate_variable(varname, varval)
        return v

    @validator("labeled_queries", each_item=True)
    def inputs_have_valid_benchkeys(cls, v: LabeledVCFInput, values) -> LabeledVCFInput:
        if "reference_sets" in values and "references" in values:
            ref_key = values["reference_sets"][v.refset].ref
            ref_benchmarks = values["references"][ref_key].benchmarks
            assert (
                v.benchmark in ref_benchmarks
            ), f"'{v.benchmark}' does not refer to a valid benchmark"
        return v

    @validator("unlabeled_queries")
    def input_keys_unique(cls, v: UnlabeledQueries, values) -> UnlabeledQueries:
        try:
            assert set(v).isdisjoint(
                set(values["labeled_queries"])
            ), "labeled and unlabeled query keys overlap overlap"
        except KeyError:
            pass
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
            tests = [t.query_key for r in v.runs.values() for t in r.test.values()]
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

    def _querykey_to_input(self, key: QueryKey) -> VCFInput:
        try:
            return self.labeled_queries[cast(LabeledQueryKey, key)]
        except KeyError:
            return self.unlabeled_queries[cast(UnlabeledQueryKey, key)]

    def refsetkey_to_refkey(self, key: RefsetKey) -> RefKey:
        return self.refsetkey_to_refset(key).ref

    def refsetkey_to_chr_indices(self, key: RefsetKey) -> Set[int]:
        f = self.refsetkey_to_refset(key).chr_filter
        return set(range(1, 25)) if len(f) == 0 else f

    def refsetkey_to_sdf_chr_filter(self, key: RefsetKey) -> Set[str]:
        indices = self.refsetkey_to_chr_indices(key)
        prefix = self.refsetkey_to_ref(key).sdf.chr_prefix
        return chr_indices_to_name(prefix, indices)

    def querykey_to_refkey(self, key: QueryKey) -> RefKey:
        rk = self.querykey_to_refsetkey(key)
        return self.refsetkey_to_refkey(rk)

    def querykey_to_refsetkey(self, key: QueryKey) -> RefsetKey:
        return self._querykey_to_input(key).refset

    def querykey_to_benchkey(self, key: LabeledQueryKey) -> BenchKey:
        return self.labeled_queries[key].benchmark

    def all_refsetkeys(self) -> Set[RefsetKey]:
        return set(
            v.refset
            for v in list(self.labeled_queries.values())
            + list(self.unlabeled_queries.values())
        )

    def all_refkeys(self) -> Set[str]:
        return set(map(self.refsetkey_to_refkey, self.all_refsetkeys()))

    def querykey_to_chr_prefix(self, key: QueryKey) -> str:
        return self._querykey_to_input(key).chr_prefix

    def querykey_to_chr_filter(self, input_key: QueryKey) -> Set[int]:
        refset_key = self.querykey_to_refsetkey(input_key)
        return self.refsetkey_to_chr_indices(refset_key)

    def querykey_to_variables(self, input_key: QueryKey) -> Dict[VarKey, str]:
        return self._querykey_to_input(input_key).variables

    def querykey_to_bench_correction(
        self,
        key: LabeledQueryKey,
    ) -> BenchmarkCorrections:
        bkey = self.labeled_queries[key].benchmark
        rkey = self.querykey_to_refsetkey(key)
        return self.refsetkey_to_ref(rkey).benchmarks[bkey].corrections

    def testkey_to_variables(
        self,
        mkey: ModelKey,
        rkey: RunKey,
        tkey: TestKey,
    ) -> Dict[VarKey, str]:
        test = self.models[mkey].runs[rkey].test[tkey]
        return {**test.variables, **self.querykey_to_variables(test.query_key)}

    def runkey_to_train_querykeys(
        self,
        mkey: ModelKey,
        rkey: RunKey,
    ) -> List[LabeledQueryKey]:
        return [t for t in self.models[mkey].runs[rkey].train]

    def runkey_to_test_querykeys(self, mkey: ModelKey, rkey: RunKey) -> List[QueryKey]:
        return [t.query_key for t in self.models[mkey].runs[rkey].test.values()]

    def testkey_to_querykey(
        self,
        mkey: ModelKey,
        rkey: RunKey,
        tkey: TestKey,
    ) -> QueryKey:
        return self.models[mkey].runs[rkey].test[tkey].query_key

    @property
    def labeled_query_resource_dir(self) -> Path:
        return self.paths.resources / "labeled_queries"

    @property
    def unlabeled_query_resource_dir(self) -> Path:
        return self.paths.resources / "unlabeled_queries"

    @property
    def ref_resource_dir(self) -> Path:
        return self.paths.resources / "reference" / all_wildcards["ref_key"]

    @property
    def bench_resource_dir(self) -> Path:
        return self.ref_resource_dir / "bench"

    @property
    def annotations_resource_dir(self) -> Path:
        return self.paths.resources / "annotations" / all_wildcards["ref_key"]

    def _result_or_log_dir(self, log: bool) -> Path:
        return self.paths.results / self.paths.log if log else self.paths.results

    def bench_dir(self, log: bool = False) -> Path:
        return (
            self._result_or_log_dir(log)
            / "bench"
            / all_wildcards["refset_key"]
            / all_wildcards["bench_key"]
        )

    def _labeled_dir(self, labeled: bool) -> Path:
        return (
            Path("labeled") / all_wildcards["l_query_key"]
            if labeled
            else Path("unlabeled") / all_wildcards["ul_query_key"]
        )

    def _query_dir(self, labeled: bool, log: bool) -> Path:
        return self._result_or_log_dir(log) / "query" / self._labeled_dir(labeled)

    def query_prepare_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "prepare"

    def query_parsed_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "parsed"

    def vcfeval_dir(self, log: bool) -> Path:
        return self._query_dir(True, log) / "prepare"

    def refset_dir(self, log: bool) -> Path:
        return self._result_or_log_dir(log) / "references" / all_wildcards["refset_key"]

    def annotated_dir(self, labeled: bool, log: bool) -> Path:
        return self._result_or_log_dir(log) / "annotated" / self._labeled_dir(labeled)

    def model_train_dir(self, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "model"
            / wildcard_format("{}-{}-{}", "model_key", "filter_key", "run_key")
        )

    def model_test_dir(self, labeled: bool, log: bool) -> Path:
        return (
            self.model_train_dir(log)
            / "test"
            / ("labeled" if labeled else "unlabeled")
            / all_wildcards["test_key"]
        )


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
        (config.querykey_to_refkey(k), v.benchmark)
        for k, v in config.labeled_queries.items()
    )
    return expand(target, zip, ref_key=rs, bench_key=bs)


class RunKeysTrain(NamedTuple):
    model_key: ModelKey
    filter_key: FilterKey
    run_key: RunKey
    labeled_query_key: LabeledQueryKey


class RunKeysTest(NamedTuple):
    model_key: ModelKey
    filter_key: FilterKey
    run_key: RunKey
    test_key: TestKey
    query_key: QueryKey


def lookup_run_sets(config: StratoMod) -> Tuple[List[RunKeysTrain], List[RunKeysTest]]:
    models = [
        ((model_key, filter_key, run_key), rest)
        for model_key, model in config.models.items()
        for filter_key in model.filter
        for run_key, rest in model.runs.items()
    ]
    train = [
        RunKeysTrain(*meta, train_key)
        for (meta, rest) in models
        for train_key in rest.train
    ]
    test = [
        RunKeysTest(*meta, test_key, test.query_key)
        for (meta, rest) in models
        for test_key, test in rest.test.items()
    ]
    return (train, test)


def test_has_bench(config: StratoMod, runs: RunKeysTest) -> bool:
    return runs.query_key in config.labeled_queries


def partition_test_set(
    config: StratoMod,
    test_set: List[RunKeysTest],
) -> Tuple[List[RunKeysTest], List[RunKeysTest]]:
    unlabeled, labeled = partition(lambda t: test_has_bench(config, t), test_set)
    return list(unlabeled), list(labeled)


def all_refset_keys(
    config: StratoMod,
    ks: Sequence[Union[RunKeysTest, RunKeysTrain]],
) -> List[RefsetKey]:
    return list(
        map(
            lambda x: config.querykey_to_refsetkey(
                x.query_key if isinstance(x, RunKeysTest) else x.labeled_query_key
            ),
            ks,
        )
    )


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
            filter_key=map(lambda x: x.filter_key.value, key_set),
            l_query_key=map(lambda x: x.labeled_query_key, key_set),
        )

    def test_targets(target: InputFiles, key_set: List[RunKeysTest]) -> InputFiles:
        return expand(
            target,
            zip,
            model_key=map(lambda x: x.model_key, key_set),
            filter_key=map(lambda x: x.filter_key.value, key_set),
            # TODO dirty hack
            l_query_key=map(lambda x: x.query_key, key_set),
            ul_query_key=map(lambda x: x.query_key, key_set),
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
            filter_key=map(lambda x: x.filter_key.value, key_set),
            run_key=map(lambda x: x.run_key, key_set),
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
        filter_key=map(lambda x: x.filter_key.value, train_set),
        run_key=map(lambda x: x.run_key, train_set),
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
