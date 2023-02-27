"""
This is the main configuration for StratoMod.

Unlike (current) snakemake pipelines, this configuration utilizes pydantic
to parse the yaml configuration into a well-defined python class. This has
several consequences/advantages:
- all python code can be statically type checked without invoking snakemake
- paths in the pipeline (which are methods in the pydantic class) can be
  resolved in the parse phase rather than at runtime.
- the user gets better errors when the configuration doesn't conform

Ideally, most errors can now be caught by linters and dry-running snakemake.

Several assumptions/conventions are maintained in order for this to work:
- all python code is segregated cleanly into conda environments, where each
  environment has mypy, pylint, etc to lint code as well as any runtime deps
- rmarkdown does not use the config at all (since R doesn't understand python
  classes)
- wildcards have consistent meaning throughout the pipeline, which allows
  matching wildcards to a given type when sent to a python script
- generally, snakemake itself is used for build dependency resolution and
  managing IO; everything else is done in raw python code which can be
  statically checked

"""
import re
import random
from pathlib import Path
from typing import (
    Generic,
    Type,
    TypeVar,
    Sequence,
    Collection,
    Callable,
    Optional,
    NamedTuple,
    NewType,
    Any,
    cast,
    Annotated,
)
from typing_extensions import Self
from more_itertools import flatten, duplicates_everseen, unzip, partition
from itertools import product
from .functional import maybe
from snakemake.io import expand, InputFiles  # type: ignore
from pydantic import BaseModel as PydanticBaseModel
from pydantic import validator, HttpUrl, PositiveFloat, NonNegativeInt, Field
from enum import Enum, auto


# various types

X = TypeVar("X")

Fraction = Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)]
NonEmptyStr = Annotated[str, Field(min_length=1)]

# newtype string wrappers for different purposes

# TODO shouldn't all these really be nonempty strings (mypy doesn't seem to
# allow this)?

RefKey = NewType("RefKey", str)  # key for a reference
RefsetKey = NewType("RefsetKey", str)  # key for reference set (reference + filter)
TestKey = NewType("TestKey", str)  # key for a model test
BenchKey = NewType("BenchKey", str)  # key for a benchmark (within a given reference)
ModelKey = NewType("ModelKey", str)  # key for a model and its parameters
RunKey = NewType("RunKey", str)  # key for a model run (eg parameters + queries)
FeatureKey = NewType("FeatureKey", str)  # key for a feature in the model

VarKey = NewType("VarKey", str)  # key for a VCF variable
VarVal = NewType("VarVal", str)  # value for a VCF variable (serialized as a string)

LabeledQueryKey = NewType("LabeledQueryKey", str)  # key for query with benchmark
UnlabeledQueryKey = NewType("UnlabeledQueryKey", str)  # key for query w/o benchmark
QueryKey = UnlabeledQueryKey | LabeledQueryKey

ChrPrefix = NewType("ChrPrefix", str)  # the "chr" (or something) prefix for chromosomes
PandasColumn = NewType("PandasColumn", str)  # the name of a pandas column

# enums to prevent runaway strings


class _ListEnum(Enum):
    "Enum which can dump all its member values easily."

    # TODO any type?
    @classmethod
    def all(cls: Type[Self]) -> list[Any]:
        return [x.value for x in cls]


class Base(_ListEnum):
    "All your bases are belong to StratoMod"
    A = "A"
    C = "C"
    G = "G"
    T = "T"


class ChrIndex(_ListEnum):
    CHR1 = auto()
    CHR2 = auto()
    CHR3 = auto()
    CHR4 = auto()
    CHR5 = auto()
    CHR6 = auto()
    CHR7 = auto()
    CHR8 = auto()
    CHR9 = auto()
    CHR10 = auto()
    CHR11 = auto()
    CHR12 = auto()
    CHR13 = auto()
    CHR14 = auto()
    CHR15 = auto()
    CHR16 = auto()
    CHR17 = auto()
    CHR18 = auto()
    CHR19 = auto()
    CHR20 = auto()
    CHR21 = auto()
    CHR22 = auto()
    CHRX = auto()
    CHRY = auto()

    @property
    def chr_name(self) -> str:
        if self.value == self.CHRY:
            return "Y"
        elif self.value == self.CHRX:
            return "X"
        return str(self.value)

    def chr_name_full(self, prefix: ChrPrefix) -> str:
        return f"{prefix}{self.chr_name}"


class FilterKey(_ListEnum):
    SNV = "SNV"
    INDEL = "INDEL"


class ErrorLabel(_ListEnum):
    FP = "fp"
    FN = "fn"


class VCFLabel(_ListEnum):
    FP = "fp"
    FN = "fn"
    TP = "tp"


class AnyLabel(_ListEnum):
    FP = "fp"
    FN = "fn"
    TP = "tp"
    TN = "tn"


class PlotType(Enum):
    BAR = "bar"
    STEP = "step"


class FeatureType(Enum):
    CONTINUOUS = "continuous"
    CATEGORICAL = "categorical"


class Binning(Enum):
    "The binning parameter in the EBM classifier"
    UNIFORM = "uniform"
    QUANTILE = "quantile"
    QUANTILE_HUMANIZED = "quantile_humanized"


class BedMergeOp(Enum):
    "Available operations for bedtools merge to perform on columns"
    MIN = "min"
    MAX = "max"
    MEAN = "mean"
    MEDIAN = "median"


class Transform(Enum):
    "Available transforms to apply to features"
    LOG = "log"
    BINARY = "binary"


# useful data bundles


class ChrFilter(NamedTuple):
    """
    Set of chromosomes and a prefix.

    Together these fully specify the desired chromosomes in a bed/ref.
    """

    prefix: ChrPrefix
    indices: set[ChrIndex]


class RunKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a model run"
    model_key: ModelKey
    filter_key: FilterKey
    run_key: RunKey


class TrainKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a train query for a run"
    run_combo: RunKeyCombo
    labeled_query_key: LabeledQueryKey


class TestKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a test query for a run"
    run_combo: RunKeyCombo
    test_key: TestKey
    query_key: QueryKey


# wildcards: all wildcards ever used in the pipeline should be in this dict


def alternate_constraint(xs: list[str]) -> str:
    return f"({'|'.join(xs)})"


# constraints should only have alphanum to avoid splitting on path slashes and
# delimiters such as "_" and "-"
KEY_CONSTR = "[a-zA-Z][a-za-z0-9]+"

_constraints: dict[str, str] = {
    # corresponds to a genome reference
    "ref_key": "[^/]+",
    # corresponds to a reference set (reference + chromosome filter + etc)
    "refset_key": "[^/]+",
    # corresponds to a testing vcf
    "test_key": "[^/]+",
    # corresponds to a query vcf file (with or without a benchmark)
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
    "label": alternate_constraint(VCFLabel.all()),
    # refers to a nucleotide base
    "base": alternate_constraint(Base.all()),
}

all_wildcards: dict[str, str] = {k: f"{{{k},{v}}}" for k, v in _constraints.items()}


def wildcard_ext(key: str, ext: str) -> str:
    return f"{all_wildcards[key]}.{ext}"


def wildcard_format(format_str: str, *keys: str) -> str:
    return format_str.format(*[all_wildcards[k] for k in keys])


def wildcard_format_ext(format_str: str, keys: str, ext: str) -> str:
    return wildcard_format(f"{format_str}.{ext}", *keys)


# snakemake config components (a pydantic model)


class _BaseModel(PydanticBaseModel):
    class Config:
        validate_all = True
        extra = "forbid"
        frozen = True


class _Range(_BaseModel):
    "Superclass for things with a lower and upper bound."
    lower: Optional[float] = None
    upper: Optional[float] = None

    @validator("upper")
    def lower_less_than_upper(
        cls,
        v: Optional[float],
        values: dict[str, Any],
    ) -> Optional[float]:
        lower = values["lower"]
        if v is not None and values["lower"] is not None:
            assert lower <= v
        return v

    def in_range(self, x: float) -> bool:
        return maybe(True, lambda y: y <= x, self.lower) and maybe(
            True, lambda y: x <= y, self.upper
        )


class Paths(_BaseModel):
    "Paths snakemake will use for disk IO"
    resources: Path
    results: Path
    log: Path


class Refset(_BaseModel):
    "A reference and a chromosome filter"
    ref: RefKey
    chr_filter: set[ChrIndex]


class CatVar(_BaseModel):
    "A categorical VCF variable"
    levels: Annotated[list[str], Field(min_items=1, unique_items=True)]


class ContVar(_Range):
    "A continuous VCF variable"
    pass


class TestDataQuery(_BaseModel):
    "Specifies a test to run within a model."
    query_key: QueryKey
    variables: dict[VarKey, VarVal]


class ModelRun(_BaseModel):
    "Specifies train and test queries to run in a model"
    train: Annotated[set[LabeledQueryKey], Field(min_items=1)]
    test: dict[TestKey, TestDataQuery]


class EBMMiscParams(_BaseModel):
    "EBM parameters that aren't specified in the classifier object"
    downsample: Optional[Fraction] = None


class EBMSplitParams(_BaseModel):
    "Parameters for the EBM train/test split"
    test_size: Fraction = 0.2
    random_state: Optional[int] = random.randrange(0, 420420)


class EBMClassifierParams(_BaseModel):
    """
    Parameters for the EBM classifier.

    The names directly correspond to those from the classifier itself.
    """

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
    random_state: int = random.randint(0, 420420)

    # for easy insertion into the ebm classifier object
    @property
    def mapping(self) -> dict[str, Any]:
        return {**self.dict(), "binning": self.binning.value}


class EBMsettings(_BaseModel):
    "All parameters for model training"
    misc_parameters: EBMMiscParams = EBMMiscParams()
    split_parameters: EBMSplitParams = EBMSplitParams()
    classifier_parameters: EBMClassifierParams = EBMClassifierParams()


class Truncation(_Range):
    "Specifies how to truncate a plot"
    pass


class Visualization(_BaseModel):
    "Specifies how to plot a feature in the final report"
    truncate: Truncation = Truncation()
    plot_type: PlotType = PlotType.STEP
    split_missing: Optional[Fraction] = None

    # convert to R-friendly dict for use in rmarkdown scripts
    @property
    def r_dict(self) -> dict[Any, Any]:
        return {**self.dict(), **{"plot_type": self.plot_type.value}}


class Feature(_BaseModel):
    "A model feature"
    feature_type: FeatureType = FeatureType.CONTINUOUS
    fill_na: Optional[float] = 0.0
    alt_name: Optional[FeatureKey] = None
    visualization: Visualization = Visualization()
    transform: Optional[Transform] = None

    # convert to R-friendly dict for use in rmarkdown scripts
    @property
    def r_dict(self) -> dict[Any, Any]:
        return {
            **self.dict(),
            **{
                "feature_type": self.feature_type.value,
                "visualization": self.visualization.r_dict,
                # pythonic fmap...
                "transform": maybe(None, lambda x: x.value, self.transform),
            },
        }


class FeaturePair(_BaseModel):
    "Two model features (for bivariate interaction terms)"
    f1: FeatureKey
    f2: FeatureKey


# NOTE need to use conlist vs set here since this will contain another
# pydantic model class, which will prevent the overall model from being
# converted to a dict (see https://github.com/pydantic/pydantic/issues/1090)
InteractionSpec_ = FeatureKey | FeaturePair
InteractionSpec = Annotated[list[InteractionSpec_], Field(unique_items=True)]


class Model(_BaseModel):
    "A fully specified model, including parameters and queries to run"
    runs: dict[RunKey, ModelRun]
    filter: set[FilterKey]
    ebm_settings: EBMsettings
    error_labels: Annotated[set[ErrorLabel], Field(min_items=1)]
    filtered_are_candidates: bool
    interactions: NonNegativeInt | InteractionSpec = 0
    features: dict[FeatureKey, Feature]

    @validator("features")
    def model_has_matching_alt_features(
        cls, fs: dict[FeatureKey, Feature]
    ) -> dict[FeatureKey, Feature]:
        for k, v in fs.items():
            alt = v.alt_name
            if alt is not None:
                prefix = assert_match("^[^_]+", k)
                alt_prefix = assert_match("^[^_]+", alt)
                assert alt_prefix == prefix, f"Alt prefix must match for {k}"
        return fs


class BedFile(_BaseModel):
    """A bed(like) file.

    'chr_prefix' must correspond to the prefix in the first column.
    """

    url: Optional[HttpUrl]
    chr_prefix: ChrPrefix


class Strats(_BaseModel):
    "Stratifications for a given reference"
    mhc: BedFile


class BenchmarkCorrections(_BaseModel):
    "Corrections to apply to a given benchmark file"
    strip_IPS: bool


class Benchmark(_BaseModel):
    "Benchmark files for a given reference"
    vcf_url: Optional[HttpUrl]
    bed_url: Optional[HttpUrl]
    chr_prefix: ChrPrefix
    corrections: BenchmarkCorrections


class Mappability(_BaseModel):
    "Low-map files for mappability features"
    low: BedFile
    high: BedFile


class Annotations(_BaseModel):
    "Bed files for feature generation"
    mappability: Mappability
    superdups: BedFile
    simreps: BedFile
    repeat_masker: BedFile


class Reference(_BaseModel):
    "A genome reference, including all associated bed/benchmark files"
    sdf: BedFile
    genome: BedFile
    strats: Strats
    annotations: Annotations
    benchmarks: dict[BenchKey, Benchmark]


FeaturePrefix = Annotated[str, Field(regex="^[A-Z]+$")]


class BedIndex(_BaseModel):
    "Metadata to track the first three columns of bed(like) dataframes"
    chr: PandasColumn
    start: PandasColumn
    end: PandasColumn

    def bed_cols_ordered(self) -> list[PandasColumn]:
        return [self.chr, self.start, self.end]

    def bed_cols_indexed(
        self, indices: tuple[int, int, int]
    ) -> dict[int, PandasColumn]:
        return dict(zip(indices, self.bed_cols_ordered()))


class _FeatureGroup(_BaseModel):
    "Superclass for feature groups (which in turn define valid feature names)"
    prefix: FeaturePrefix

    def fmt_feature(self, rest: str) -> FeatureKey:
        return FeatureKey(f"{self.prefix}_{rest}")


class VCFColumns(_BaseModel):
    "Columns for a vcf file"
    qual: NonEmptyStr
    filter: NonEmptyStr
    info: NonEmptyStr
    gt: NonEmptyStr
    gq: NonEmptyStr
    dp: NonEmptyStr
    vaf: NonEmptyStr
    len: NonEmptyStr


class VCFGroup(_FeatureGroup):
    "Feature and column names for VCF files"
    columns: VCFColumns

    def fmt_name(self: Self, f: Callable[[VCFColumns], str]) -> FeatureKey:
        return self.fmt_feature(f(self.columns))

    @property
    def str_feature_names(self) -> list[FeatureKey]:
        return [
            self.fmt_name(x)
            for x in [
                lambda x: x.filter,
                lambda x: x.info,
                lambda x: x.gt,
                lambda x: x.gq,
            ]
        ]

    def feature_names(self) -> list[FeatureKey]:
        return self.str_feature_names + [
            self.fmt_name(x)
            for x in [
                lambda x: x.qual,
                lambda x: x.dp,
                lambda x: x.vaf,
                lambda x: x.len,
            ]
        ]


class MapSuffixes(_BaseModel):
    "Suffixes corresponding to low-map regions"
    low: NonEmptyStr
    high: NonEmptyStr


class MapGroup(_FeatureGroup):
    "Feature and column names for low-map dataframes"
    suffixes: MapSuffixes

    @property
    def low(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.low)

    @property
    def high(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.high)

    def feature_names(self) -> list[FeatureKey]:
        return [self.low, self.high]


class HomopolySuffixes(_BaseModel):
    "Suffixes corresponding to homopolymer regions"
    len: NonEmptyStr
    imp_frac: NonEmptyStr


class HomopolyGroup(_FeatureGroup):
    "Feature and column names for homopolymer dataframes"
    bases: set[Base]
    suffixes: HomopolySuffixes

    def fmt_name(self, b: Base, f: Callable[[HomopolySuffixes], str]) -> FeatureKey:
        return self.fmt_feature(f"{b.value}_{f(self.suffixes)}")

    def feature_names(self) -> list[FeatureKey]:
        return [
            self.fmt_name(b, f)
            for f, b in product([lambda x: x.len, lambda x: x.imp_frac], self.bases)
        ]


class RMSKGroup(_FeatureGroup):
    "Feature and column names for repeat masker dataframes"
    classes: dict[str, set[NonEmptyStr]] = {
        k: set() for k in ["SINE", "LINE", "LTR", "Satellite"]
    }

    def fmt_name(
        self,
        grp: str,
        fam: Optional[str],
    ) -> PandasColumn:
        assert grp in self.classes, f"{grp} not a valid RMSK class"
        rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
        return PandasColumn(self.fmt_feature(f"{rest}_length"))

    def feature_names(self) -> list[FeatureKey]:
        def fmt(grp: str, fam: Optional[str]) -> FeatureKey:
            return FeatureKey(self.fmt_name(grp, fam))

        return [
            *[fmt(c, None) for c in self.classes],
            *[fmt(c, f) for c, fs in self.classes.items() for f in fs],
        ]


class MergedFeatureGroup(_FeatureGroup, Generic[X]):
    "Superclass for feature group which supports bedtools merge"
    operations: set[BedMergeOp]
    columns: X

    def fmt_col(self, f: Callable[[X], str]) -> PandasColumn:
        return PandasColumn(self.fmt_feature(f(self.columns)))

    def fmt_count_feature(self) -> FeatureKey:
        return self.fmt_feature("count")

    def fmt_merged_feature(self, middle: str, op: BedMergeOp) -> FeatureKey:
        return FeatureKey(f"{middle}_{op.value}")

    def merged_feature_col_names(
        self,
        fs: list[Callable[[X], str]],
    ) -> list[FeatureKey]:
        return self.merged_feature_names([self.fmt_col(f) for f in fs])

    def merged_feature_names(self, names: list[PandasColumn]) -> list[FeatureKey]:
        return [
            *[
                self.fmt_merged_feature(n, o)
                for n, o in product(names, self.operations)
            ],
            self.fmt_count_feature(),
        ]


class SegDupsColumns(_BaseModel):
    "Columns corresponding to the superdups files"
    alignL: NonEmptyStr
    fracMatchIndel: NonEmptyStr


class SegDupsGroup(MergedFeatureGroup[SegDupsColumns]):
    "Feature and column names for segdups dataframes"
    columns: SegDupsColumns

    def feature_names(self) -> list[FeatureKey]:
        return self.merged_feature_col_names(
            [lambda x: x.alignL, lambda x: x.fracMatchIndel]
        )


class TandemRepeatColumns(_BaseModel):
    "Columns corresponding to the simple_repeats files"
    period: NonEmptyStr
    copyNum: NonEmptyStr
    perMatch: NonEmptyStr
    perIndel: NonEmptyStr
    score: NonEmptyStr


class TandemRepeatGroup(MergedFeatureGroup[TandemRepeatColumns]):
    "Feature and column names for segdups dataframes"
    columns: TandemRepeatColumns

    def _base_name(self, base: str) -> PandasColumn:
        return PandasColumn(self.fmt_feature(f"percent_{base}"))

    @property
    def AT_name(self) -> PandasColumn:
        return self._base_name("AT")

    @property
    def GC_name(self) -> PandasColumn:
        return self._base_name("GC")

    @property
    def AG_name(self) -> PandasColumn:
        return self._base_name("AG")

    @property
    def CT_name(self) -> PandasColumn:
        return self._base_name("CT")

    @property
    def length_name(self) -> PandasColumn:
        return PandasColumn(self.fmt_feature("length"))

    def fmt_base_col(self, b: Base) -> PandasColumn:
        return self._base_name(b.value)

    def feature_names(self) -> list[FeatureKey]:
        single_base = [self.fmt_base_col(b) for b in Base]
        # lombardo quit to become a bioinformatician...
        double_base = [self.AT_name, self.GC_name, self.AG_name, self.CT_name]
        noncols = self.merged_feature_names(
            single_base + double_base + [self.length_name]
        )
        cols = self.merged_feature_col_names(
            [
                lambda x: x.period,
                lambda x: x.copyNum,
                lambda x: x.perMatch,
                lambda x: x.perIndel,
                lambda x: x.score,
            ]
        )
        return noncols + cols


class VariableGroup(_FeatureGroup):
    "Feature and column names for vcf manually-assigned variable dataframes"
    continuous: dict[VarKey, ContVar]
    categorical: dict[VarKey, CatVar]

    # not a pydantic validator; requires lots of data from parent classes
    def validate_variable(self, varname: VarKey, varval: VarVal) -> None:
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

    def all_keys(self) -> list[VarKey]:
        return list(self.continuous) + list(self.categorical)

    def feature_names(self) -> list[FeatureKey]:
        return [self.fmt_feature(x) for x in self.all_keys()]

    # ASSUME we don't need to check the incoming varkeys or varvals since they
    # will be validated on model creation
    def _parse_var(self, k: VarKey, v: VarVal) -> float:
        if k in self.categorical:
            return self.categorical[k].levels.index(v)
        elif k in self.continuous:
            return float(v)
        else:
            assert False, f"{k} not a valid variable key: this should not happen"

    def parse_vars(self, kvs: dict[VarKey, VarVal]) -> dict[FeatureKey, float]:
        return {self.fmt_feature(k): self._parse_var(k, v) for k, v in kvs.items()}


class FormatFields(_BaseModel):
    "Specifies the names of fields in the FORMAT column of a VCF file"
    vaf: Optional[str] = "VAF"
    dp: Optional[str] = "DP"
    gt: Optional[str] = "GT"
    gq: Optional[str] = "GQ"

    def vcf_fields(self, vcf: VCFGroup) -> dict[PandasColumn, Optional[str]]:
        return {
            PandasColumn(vcf.fmt_name(f)): v
            for f, v in [
                (lambda x: x.vaf, self.vaf),
                (lambda x: x.dp, self.dp),
                (lambda x: x.gt, self.gt),
                (lambda x: x.gq, self.gq),
            ]
        }


class UnlabeledVCFQuery(_BaseModel):
    "A vcf to be used as the query for a model without labels."
    refset: RefsetKey
    chr_prefix: ChrPrefix
    url: Optional[HttpUrl]
    variables: dict[VarKey, VarVal]
    format_fields: FormatFields = FormatFields()
    max_ref: Annotated[int, Field(ge=0)] = 50
    max_alt: Annotated[int, Field(ge=0)] = 50


class LabeledVCFQuery(UnlabeledVCFQuery):
    "A vcf to be used as the query for a model with labels (requires benchmark)."
    benchmark: BenchKey


VCFQuery = UnlabeledVCFQuery | LabeledVCFQuery


class FeatureNames(_BaseModel):
    "Defines valid feature names to be specified in models"
    label: NonEmptyStr
    raw_index: NonEmptyStr
    bed_index: BedIndex
    vcf: VCFGroup
    mappability: MapGroup
    homopolymers: HomopolyGroup
    repeat_masker: RMSKGroup
    segdups: SegDupsGroup
    tandem_repeats: TandemRepeatGroup
    variables: VariableGroup

    @property
    def label_name(self: Self) -> PandasColumn:
        return PandasColumn(self.label)

    def all_index_cols(self) -> list[PandasColumn]:
        return [PandasColumn(self.raw_index), *self.bed_index.bed_cols_ordered()]

    def all_feature_names(self) -> set[FeatureKey]:
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

    @property
    def non_summary_cols(self) -> list[PandasColumn]:
        return self.all_index_cols() + [
            PandasColumn(x) for x in self.vcf.str_feature_names
        ]


class Tools(_BaseModel):
    "Various tools that need to be downloaded for the pipeline to work"
    repseq: HttpUrl


LabeledQueries = dict[LabeledQueryKey, LabeledVCFQuery]
UnlabeledQueries = dict[UnlabeledQueryKey, UnlabeledVCFQuery]


def assert_empty(xs: Collection[X], msg: str) -> None:
    assert len(xs) == 0, f"{msg}: {xs}"


def assert_no_dups(xs: Collection[X], msg: str) -> None:
    assert_empty(set(duplicates_everseen(xs)), msg)


def assert_match(pat: str, s: str) -> str:
    res = re.match(pat, s)
    assert res is not None, f"match failed for pattern '{pat}' and query '{s}'"
    return res[0]


def assert_subset(xs: set[X], ys: set[X]) -> None:
    assert xs <= ys, f"not a subset - extra members: {xs - ys}"


def flatten_features(fs: dict[FeatureKey, Feature]) -> list[FeatureKey]:
    return [k if v.alt_name is None else v.alt_name for k, v in fs.items()]


class StratoMod(_BaseModel):
    "Root config for the stratomod pipeline."

    paths: Paths
    tools: Tools
    references: dict[RefKey, Reference]
    feature_names: FeatureNames
    reference_sets: dict[RefsetKey, Refset]
    labeled_queries: LabeledQueries
    unlabeled_queries: UnlabeledQueries
    models: dict[ModelKey, Model]

    @validator("reference_sets", each_item=True)
    def refsets_have_valid_refkeys(
        cls: Type[Self],
        v: Refset,
        values: dict[str, Any],
    ) -> Refset:
        try:
            assert (
                v.ref in values["references"]
            ), f"'{v.ref}' does not refer to a valid reference"
        except KeyError:
            pass
        return v

    @validator("labeled_queries", "unlabeled_queries", each_item=True)
    def inputs_have_valid_refsetkeys(
        cls: Type[Self],
        v: VCFQuery,
        values: dict[str, Any],
    ) -> VCFQuery:
        try:
            assert (
                v.refset in values["reference_sets"]
            ), f"'{v.refset}' does not refer to a valid reference set"
        except KeyError:
            pass
        return v

    @validator("labeled_queries", "unlabeled_queries", each_item=True)
    def inputs_have_valid_variables(
        cls: Type[Self],
        v: VCFQuery,
        values: dict[str, Any],
    ) -> VCFQuery:
        try:
            var_root = cast(FeatureNames, values["feature_names"]).variables
        except KeyError:
            pass
        else:
            for varname, varval in v.variables.items():
                var_root.validate_variable(varname, varval)
        return v

    @validator("labeled_queries", each_item=True)
    def inputs_have_valid_benchkeys(
        cls: Type[Self],
        v: LabeledVCFQuery,
        values: dict[str, Any],
    ) -> LabeledVCFQuery:
        try:
            refsets = values["reference_sets"]
            refs = values["references"]
        except KeyError:
            pass
        else:
            ref_key = refsets[v.refset].ref
            ref_benchmarks = refs[ref_key].benchmarks
            assert (
                v.benchmark in ref_benchmarks
            ), f"'{v.benchmark}' does not refer to a valid benchmark"
        return v

    @validator("unlabeled_queries")
    def input_keys_unique(
        cls: Type[Self],
        v: UnlabeledQueries,
        values: dict[str, Any],
    ) -> UnlabeledQueries:
        try:
            assert set(v).isdisjoint(
                set(values["labeled_queries"])
            ), "labeled and unlabeled query keys overlap"
        except KeyError:
            pass
        return v

    @validator("models", each_item=True)
    def models_have_valid_features(
        cls: Type[Self],
        v: Model,
        values: dict[str, Any],
    ) -> Model:
        try:
            features = values["feature_names"].all_feature_names()
        except KeyError:
            pass
        else:
            assert_subset(set(v.features), features)
        return v

    @validator("models", each_item=True)
    def models_have_valid_features_alt(cls: Type[Self], v: Model) -> Model:
        # TODO dry?
        assert_no_dups(flatten_features(v.features), "Duplicated features")
        return v

    @validator("models", each_item=True)
    def models_have_valid_interactions(
        cls: Type[Self],
        v: Model,
        values: dict[str, Any],
    ) -> Model:
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
    def models_have_valid_runs_train(
        cls: Type[Self],
        v: Model,
        values: dict[str, Any],
    ) -> Model:
        try:
            train = [t for r in v.runs.values() for t in r.train]
            assert_subset(set(train), set(values["inputs"]))
        except KeyError:
            pass
        return v

    @validator("models", each_item=True)
    def models_have_valid_runs_test(
        cls: Type[Self],
        v: Model,
        values: dict[str, Any],
    ) -> Model:
        try:
            tests = [t.query_key for r in v.runs.values() for t in r.test.values()]
            assert_subset(set(tests), set(values["inputs"]))
        except KeyError:
            pass
        return v

    # TODO might be nice to alert user when they try and test a variable that's
    # not included in the train set, since this will likely screw with things
    # in weird ways

    @validator("models", each_item=True)
    def models_have_valid_runs_test_variables(
        cls: Type[Self],
        v: Model,
        values: dict[str, Any],
    ) -> Model:
        try:
            var_root = cast(FeatureNames, values["feature_names"]).variables
        except KeyError:
            pass
        else:
            varpairs = [
                (varname, varval)
                for r in v.runs.values()
                for t in r.test.values()
                for varname, varval in t.variables.items()
            ]
            assert set([p[0] for p in varpairs]) == set(
                var_root.all_keys()
            ), "missing variables"
            # assert_subset(set([p[0] for p in varpairs]), set(var_root.all_keys()))
            for varname, varval in varpairs:
                var_root.validate_variable(varname, varval)
        return v

    # various mapping functions

    def refsetkey_to_ref(self, key: RefsetKey) -> Reference:
        return self.references[self.refsetkey_to_refkey(key)]

    def refsetkey_to_refset(self, key: RefsetKey) -> Refset:
        return self.reference_sets[key]

    def _querykey_to_input(self, key: QueryKey) -> VCFQuery:
        try:
            return self.labeled_queries[cast(LabeledQueryKey, key)]
        except KeyError:
            return self.unlabeled_queries[cast(UnlabeledQueryKey, key)]

    def refsetkey_to_refkey(self, key: RefsetKey) -> RefKey:
        return self.refsetkey_to_refset(key).ref

    def refsetkey_to_chr_indices(self, key: RefsetKey) -> set[ChrIndex]:
        f = self.refsetkey_to_refset(key).chr_filter
        return set(x for x in ChrIndex) if len(f) == 0 else f

    def refsetkey_to_chr_filter(
        self,
        get_prefix: Callable[[Reference], ChrPrefix],
        key: RefsetKey,
    ) -> ChrFilter:
        indices = self.refsetkey_to_chr_indices(key)
        prefix = get_prefix(self.refsetkey_to_ref(key))
        return ChrFilter(prefix, indices)

    def refsetkey_to_sdf_chr_filter(self, key: RefsetKey) -> set[str]:
        prefix, indices = self.refsetkey_to_chr_filter(lambda r: r.sdf.chr_prefix, key)
        return set(i.chr_name_full(prefix) for i in indices)

    def benchkey_to_chr_prefix(self, rkey: RefsetKey, bkey: BenchKey) -> str:
        return self.refsetkey_to_ref(rkey).benchmarks[bkey].chr_prefix

    def refkey_to_annotations(self, key: RefKey) -> Annotations:
        return self.references[key].annotations

    def querykey_to_refkey(self, key: QueryKey) -> RefKey:
        rk = self.querykey_to_refsetkey(key)
        return self.refsetkey_to_refkey(rk)

    def querykey_to_refsetkey(self, key: QueryKey) -> RefsetKey:
        return self._querykey_to_input(key).refset

    def querykey_to_benchkey(self, key: LabeledQueryKey) -> BenchKey:
        return self.labeled_queries[key].benchmark

    def querykey_to_chr_prefix(self, key: QueryKey) -> str:
        return self._querykey_to_input(key).chr_prefix

    def querykey_to_variables(self, input_key: QueryKey) -> dict[FeatureKey, float]:
        vs = self._querykey_to_input(input_key).variables
        return self.feature_names.variables.parse_vars(vs)

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
    ) -> dict[FeatureKey, float]:
        test = self.models[mkey].runs[rkey].test[tkey]
        qs = self.querykey_to_variables(test.query_key)
        return {**qs, **self.feature_names.variables.parse_vars(test.variables)}

    def runkey_to_train_querykeys(
        self,
        mkey: ModelKey,
        rkey: RunKey,
    ) -> list[LabeledQueryKey]:
        return [t for t in self.models[mkey].runs[rkey].train]

    def runkey_to_test_querykeys(self, mkey: ModelKey, rkey: RunKey) -> list[QueryKey]:
        return [t.query_key for t in self.models[mkey].runs[rkey].test.values()]

    def testkey_to_querykey(
        self,
        mkey: ModelKey,
        rkey: RunKey,
        tkey: TestKey,
    ) -> QueryKey:
        return self.models[mkey].runs[rkey].test[tkey].query_key

    # path expansion

    def _workflow_path(self, components: list[str]) -> Path:
        p = Path(*components).resolve()
        assert p.exists(), f"{p} does not exist"
        return p

    def env_file(self, envname: str) -> Path:
        return self._workflow_path(["workflow/envs", f"{envname}.yml"])

    def _scripts_dir(self, rest: list[str]) -> Path:
        return self._workflow_path(["workflow/scripts", *rest])

    def python_script(self, basename: str) -> Path:
        return self._scripts_dir(["python", basename])

    def rmd_script(self, basename: str) -> Path:
        return self._scripts_dir(["rmarkdown", basename])

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

    def annotation_resource_dir(self, which: str) -> Path:
        return self.paths.resources / "annotations" / all_wildcards["ref_key"] / which

    @property
    def tool_resource_dir(self) -> Path:
        return self.paths.resources / "tools"

    def tool_dir(self, log: bool) -> Path:
        return self._result_or_log_dir(log) / "tools"

    def _result_or_log_dir(self, log: bool) -> Path:
        return self.paths.results / self.paths.log if log else self.paths.results

    def bench_dir(self, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "bench"
            / all_wildcards["refset_key"]
            / all_wildcards["bench_key"]
        )

    def annotation_dir(self, which: str, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "annotations"
            / all_wildcards["refset_key"]
            / which
        )

    def _labeled_dir(self, labeled: bool) -> Path:
        return (
            Path("labeled") / all_wildcards["l_query_key"]
            if labeled
            else Path("unlabeled") / all_wildcards["ul_query_key"]
        )

    def _query_dir(self, labeled: bool, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "query"
            / all_wildcards["refset_key"]
            / self._labeled_dir(labeled)
        )

    def query_prepare_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "prepare"

    def query_parsed_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "parsed"

    def vcfeval_dir(self, log: bool) -> Path:
        return self._query_dir(True, log) / "vcfeval"

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

    # wildcard expansion

    def querykey_is_labeled(self, k: QueryKey) -> bool:
        return k in self.labeled_queries

    @property
    def _all_runs(self) -> list[ModelRun]:
        return [r for m in self.models.values() for r in m.runs.values()]

    @property
    def all_labeled_querykeys(self) -> set[LabeledQueryKey]:
        runs = self._all_runs
        train_keys = [k for r in runs for k in r.train]
        test_keys = [
            LabeledQueryKey(k)
            for r in runs
            for t in r.test.values()
            if self.querykey_is_labeled(k := t.query_key)
        ]
        return set(train_keys + test_keys)

    @property
    def all_unlabeled_querykeys(self) -> set[UnlabeledQueryKey]:
        return set(
            [
                UnlabeledQueryKey(k)
                for r in self._all_runs
                for t in r.test.values()
                if not self.querykey_is_labeled(k := t.query_key)
            ]
        )

    @property
    def all_refsetkeys(self) -> set[RefsetKey]:
        return set(
            list(map(self.querykey_to_refsetkey, self.all_labeled_querykeys))
            + list(map(self.querykey_to_refsetkey, self.all_unlabeled_querykeys))
        )

    @property
    def all_refkeys(self) -> set[str]:
        return set(map(self.refsetkey_to_refkey, self.all_refsetkeys))

    @property
    def _all_model_combos(self) -> tuple[list[TrainKeyCombo], list[TestKeyCombo]]:
        models = [
            (RunKeyCombo(model_key, filter_key, run_key), rest)
            for model_key, model in self.models.items()
            for filter_key in model.filter
            for run_key, rest in model.runs.items()
        ]
        train = [
            TrainKeyCombo(r, train_key)
            for (r, rest) in models
            for train_key in rest.train
        ]
        test = [
            TestKeyCombo(r, test_key, test.query_key)
            for (r, rest) in models
            for test_key, test in rest.test.items()
        ]
        return (train, test)

    def bench_targets(self: Self, target: InputFiles) -> InputFiles:
        rs, bs = unzip(
            (self.querykey_to_refkey(k), v.benchmark)
            for k, v in self.labeled_queries.items()
        )
        return expand(target, zip, ref_key=rs, bench_key=bs)

    def partition_test_set(
        self,
        test_set: list[TestKeyCombo],
    ) -> tuple[list[TestKeyCombo], list[TestKeyCombo]]:
        unlabeled, labeled = partition(
            lambda t: self.querykey_is_labeled(t.query_key), test_set
        )
        return list(unlabeled), list(labeled)

    def summary_targets(
        self,
        labeled_target: InputFiles,
        unlabeled_target: InputFiles,
    ) -> InputFiles:
        def labeled_targets(
            target: InputFiles, key_set: list[TrainKeyCombo]
        ) -> InputFiles:
            return expand(
                target,
                zip,
                model_key=map(lambda x: x.run_combo.model_key, key_set),
                filter_key=map(lambda x: x.run_combo.filter_key.value, key_set),
                l_query_key=map(lambda x: x.labeled_query_key, key_set),
            )

        def test_targets(
            target: InputFiles,
            key_set: list[TestKeyCombo],
            labeled: bool,
        ) -> InputFiles:
            qkey = "l_query_key" if labeled else "ul_query_key"
            return expand(
                target,
                zip,
                model_key=map(lambda x: x.run_combo.model_key, key_set),
                filter_key=map(lambda x: x.run_combo.filter_key.value, key_set),
                **{qkey: map(lambda x: x.query_key, key_set)},
            )

        train_set, test_set = self._all_model_combos
        unlabeled_test_set, labeled_test_set = self.partition_test_set(test_set)

        return (
            labeled_targets(labeled_target, train_set)
            + test_targets(labeled_target, labeled_test_set, True)
            + test_targets(unlabeled_target, unlabeled_test_set, False)
        )

    def model_targets(
        self,
        train_target: InputFiles,
        labeled_test_target: InputFiles,
        unlabeled_test_target: InputFiles,
    ) -> InputFiles:
        def all_refset_keys(
            ks: Sequence[TestKeyCombo | TrainKeyCombo],
        ) -> list[RefsetKey]:
            return list(
                map(
                    lambda x: self.querykey_to_refsetkey(
                        x.query_key
                        if isinstance(x, TestKeyCombo)
                        else x.labeled_query_key
                    ),
                    ks,
                )
            )

        def test_targets(
            path: InputFiles,
            key_set: list[TestKeyCombo],
            labeled: bool,
        ) -> InputFiles:
            qkey = "l_query_key" if labeled else "ul_query_key"
            return expand(
                path,
                zip,
                model_key=map(lambda x: x.run_combo.model_key, key_set),
                filter_key=map(lambda x: x.run_combo.filter_key.value, key_set),
                run_key=map(lambda x: x.run_combo.run_key, key_set),
                test_key=map(lambda x: x.test_key, key_set),
                refset_key=all_refset_keys(train_set),
                **{qkey: map(lambda x: x.query_key, key_set)},
            )

        train_set, test_set = self._all_model_combos
        unlabeled_test_set, labeled_test_set = self.partition_test_set(test_set)
        train = expand(
            train_target,
            zip,
            model_key=map(lambda x: x.run_combo.model_key, train_set),
            filter_key=map(lambda x: x.run_combo.filter_key.value, train_set),
            run_key=map(lambda x: x.run_combo.run_key, train_set),
        )

        # TODO these should eventually point to the test summary htmls
        labeled_test = test_targets(labeled_test_target, labeled_test_set, True)
        unlabeled_test = test_targets(unlabeled_test_target, unlabeled_test_set, False)

        return train + labeled_test + unlabeled_test

    # NOTE: dirty hack to get this to work with rmarkdown; since I don't use
    # the config in rmarkdown, just return a blank dict. Returning a real dict
    # will work assuming all the types in the underlying structure are hashable
    # but not all will convert to R types (eg enum, which is silly)
    def items(self: Self) -> Any:
        return {}.items()


# other random functions


def attempt_mem_gb(mem_gb: int) -> Callable[[dict[str, str], int], int]:
    # double initial memory on each attempt
    return lambda _, attempt: cast(int, mem_gb * 1000 * 2 ** (attempt - 1))
