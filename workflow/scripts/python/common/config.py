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
from pydantic import validator, HttpUrl, PositiveFloat, NonNegativeInt, Field, FilePath
from enum import Enum, unique
from collections import ChainMap


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
FeatureKey = NewType("FeatureKey", str)  # key for a feature in the model

VarKey = NewType("VarKey", str)  # key for a VCF variable
VarVal = NewType("VarVal", str)  # value for a VCF variable (serialized as a string)

LabeledQueryKey = NewType("LabeledQueryKey", str)  # key for query with benchmark
UnlabeledQueryKey = NewType("UnlabeledQueryKey", str)  # key for query w/o benchmark
QueryKey = UnlabeledQueryKey | LabeledQueryKey

ChrPrefix = NewType("ChrPrefix", str)  # the "chr" (or something) prefix for chromosomes
PandasColumn = NewType("PandasColumn", str)  # the name of a pandas column
FeatureDesc = NewType("FeatureDesc", str)  # a description for a feature

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


@unique
class ChrIndex(Enum):
    _ignore_ = "ChrIndex i"
    ChrIndex = vars()
    for i in range(1, 23):
        ChrIndex[f"CHR{i}"] = i
    CHRX = 23
    CHRY = 24

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, prefix: str) -> str:
        return f"{prefix}{self.chr_name}"

    def __lt__(self, other: Self) -> bool:
        return cast(bool, self.value < other.value)


class VartypeKey(_ListEnum):
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


# constants

BED_CHROM = "chrom"
BED_START = "chromStart"
BED_END = "chromEnd"
VAR_IDX = "variant_index"

BED_COLS = [BED_CHROM, BED_START, BED_END]
IDX_COLS = [VAR_IDX, *BED_COLS]


def bed_cols_indexed(indices: tuple[int, int, int]) -> dict[int, str]:
    return dict(zip(indices, BED_COLS))


# useful data bundles


class ChrFilter(NamedTuple):
    """
    Set of chromosomes and a prefix.

    Together these fully specify the desired chromosomes in a bed/ref.

    The indices are assumed to be ordered.
    """

    prefix: ChrPrefix
    indices: list[ChrIndex]


class ModelKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a model run"
    model_key: ModelKey
    vartype_key: VartypeKey


class TrainKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a train query for a run"
    run_combo: ModelKeyCombo
    labeled_query_key: LabeledQueryKey


class TestKeyCombo(NamedTuple):
    "The keys needed to uniquely specify a test query for a run"
    run_combo: ModelKeyCombo
    test_key: TestKey
    query_key: QueryKey


# wildcards: all wildcards ever used in the pipeline should be in this dict


def alternate_constraint(xs: list[str]) -> str:
    return f"({'|'.join(xs)})"


# constraints should only have alphanum or underscore to avoid splitting on path
# delimitors ("/"), file extensions ("."), or filename delimiters ("-" or
# similar)
KEY_CONSTR = "[A-Za-z][A-Za-z0-9_]+"

_constraints: dict[str, str] = {
    # corresponds to a genome reference
    "ref_key": KEY_CONSTR,
    # corresponds to a reference set (reference + chromosome filter + etc)
    "refset_key": KEY_CONSTR,
    # corresponds to a testing vcf
    "test_key": KEY_CONSTR,
    # corresponds to a query vcf file (with or without a benchmark)
    "query_key": KEY_CONSTR,
    "ul_query_key": KEY_CONSTR,
    "l_query_key": KEY_CONSTR,
    # refers to a benchmark vcf (within the context of a given reference)
    "bench_key": KEY_CONSTR,
    # refers to an EBM model and its parameters and inputs
    "model_key": KEY_CONSTR,
    # refers to the variant type (SNP or INDEL, for now)
    "vartype_key": alternate_constraint(VartypeKey.all()),
    # refers to a variant benchmarking label (tp, fp, etc)
    "label": alternate_constraint(VCFLabel.all()),
    # refers to a nucleotide base
    "base": alternate_constraint(Base.all()),
}


def get_wildcard(x: str) -> str:
    assert x in _constraints, f"invalid wildcard: {x}"
    return f"{{{x}}}"


def wildcard_ext(key: str, ext: str) -> str:
    return f"{get_wildcard(key)}.{ext}"


def wildcard_format(format_str: str, *keys: str) -> str:
    return format_str.format(*[get_wildcard(k) for k in keys])


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
    lower: float | None = None
    upper: float | None = None

    @validator("upper")
    def lower_less_than_upper(
        cls,
        v: float | None,
        values: dict[str, Any],
    ) -> float | None:
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
    description: FeatureDesc | None = None


class ContVar(_Range):
    "A continuous VCF variable"
    description: FeatureDesc | None = None


class TestDataQuery(_BaseModel):
    "Specifies a test to run within a model."
    query_key: QueryKey
    variables: dict[VarKey, VarVal] = {}


class EBMMiscParams(_BaseModel):
    "EBM parameters that aren't specified in the classifier object"
    downsample: Fraction | None = None


class EBMSplitParams(_BaseModel):
    "Parameters for the EBM train/test split"
    test_size: Fraction = 0.2
    random_state: int | None = random.randrange(0, 420420)


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
    split_missing: Fraction | None = None

    # convert to R-friendly dict for use in rmarkdown scripts
    @property
    def r_dict(self) -> dict[Any, Any]:
        return {**self.dict(), **{"plot_type": self.plot_type.value}}


class Feature(_BaseModel):
    "A model feature"
    feature_type: FeatureType = FeatureType.CONTINUOUS
    fill_na: float | None = 0.0
    alt_name: FeatureKey | None = None
    visualization: Visualization = Visualization()
    transform: Transform | None = None

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
    train: Annotated[set[LabeledQueryKey], Field(min_items=1)]
    test: dict[TestKey, TestDataQuery]
    vartypes: set[VartypeKey]
    ebm_settings: EBMsettings = EBMsettings()
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


class HashedSrc(_BaseModel):
    md5: str | None = None


class LocalSrc(HashedSrc):
    filepath: FilePath


class HTTPSrc(HashedSrc):
    url: HttpUrl


FileSrc = LocalSrc | HTTPSrc


class BedColumns(_BaseModel):
    """Denotes coordinate columns in a bed file (0-indexed)."""

    chrom: NonNegativeInt = 0
    start: NonNegativeInt = 1
    end: NonNegativeInt = 2

    @validator("start")
    def start_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert values["chr"] != v, "Bed columns must be different"
        except KeyError:
            pass
        return v

    @validator("end")
    def end_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert (
                values["chr"] != v and values["start"] != v
            ), "Bed columns must be different"
        except KeyError:
            pass
        return v

    def assert_different(self, x: int) -> None:
        assert (
            self.chrom != x and self.start != x and self.end != x
        ), "Column must be different index"

    @property
    def typed(self) -> dict[int, Type[int | str]]:
        return {self.chrom: str, self.start: int, self.end: int}

    @property
    def indexed(self) -> dict[int, str]:
        return bed_cols_indexed((self.chrom, self.start, self.end))


class VCFFile(_BaseModel):
    """A VCF file.

    'chr_prefix' must correspond to the prefix in the first column.
    """

    src: FileSrc
    chr_prefix: ChrPrefix = ChrPrefix("chr")


class BedFileParams(_BaseModel):
    """Parameters decribing how to parse a bed-like file.

    Members:
    chr_prefix - the prefix on the chromosomes
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    chr_prefix: ChrPrefix = ChrPrefix("chr")
    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"


class BedFile(_BaseModel):
    """A bed(like) file."""

    src: FileSrc
    params: BedFileParams = BedFileParams()


class TRFColumns(_BaseModel):
    period: NonNegativeInt = 5
    copy_num: NonNegativeInt = 6
    per_match: NonNegativeInt = 8
    per_indel: NonNegativeInt = 9
    score: NonNegativeInt = 10
    per_A: NonNegativeInt = 11
    per_C: NonNegativeInt = 12
    per_G: NonNegativeInt = 13
    per_T: NonNegativeInt = 14


class TandemRepeatsFile(BedFile):
    """A tandem-repeat-finder output file (or similar)"""

    other_cols: TRFColumns = TRFColumns()


class SegdupsColumns(_BaseModel):
    align_L: NonNegativeInt = 18
    frac_match_indel: NonNegativeInt = 27


class SegdupsFile(BedFile):
    """A superdups-like file"""

    other_cols: SegdupsColumns = SegdupsColumns()


class RMSKColumns(_BaseModel):
    cls: NonNegativeInt = 11
    family: NonNegativeInt = 12


class RMSKFile(BedFile):
    """A repeat-masker-like output file"""

    other_cols: RMSKColumns = RMSKColumns()
    class_families: dict[str, list[str]]


class RefFile(_BaseModel):
    src: FileSrc
    is_fasta: bool
    chr_prefix: ChrPrefix = ChrPrefix("chr")


class BedRegion(_BaseModel):
    chrom: ChrIndex
    start: NonNegativeInt
    end: NonNegativeInt

    @validator("end")
    def positive_region(cls, v: int, values: dict[Any, Any]) -> int:
        try:
            start = cast(BedRegion, values["feature_definitions"]).start
        except KeyError:
            pass
        else:
            assert v > start, "End must be greater than start"
        return v

    def fmt(self) -> str:
        return "\t".join(map(str, [self.chrom.value, self.start, self.start]))

    def __lt__(self, other: Self) -> bool:
        return (
            self.chrom,
            self.start,
            self.end,
        ) < (
            other.chrom,
            other.start,
            other.end,
        )


class Strats(_BaseModel):
    "Stratifications for a given reference"
    mhc: list[BedRegion]


class BenchmarkCorrections(_BaseModel):
    "Corrections to apply to a given benchmark file"
    strip_IPS: bool


class Benchmark(_BaseModel):
    "Benchmark files for a given reference"
    vcf: VCFFile
    bed: BedFile
    corrections: BenchmarkCorrections


class Mappability(_BaseModel):
    "Low-map files for mappability features"
    low: BedFile
    high: BedFile


class FeatureData(_BaseModel):
    "Bed files for feature generation"
    mappability: Mappability
    segdups: SegdupsFile
    tandem_repeats: TandemRepeatsFile
    repeat_masker: RMSKFile


class Reference(_BaseModel):
    "A genome reference, including all associated bed/benchmark files"
    sdf: RefFile
    strats: Strats
    feature_data: FeatureData
    benchmarks: dict[BenchKey, Benchmark]


FeaturePrefix = Annotated[str, Field(regex="^[A-Z]+$")]


class _FeatureGroup(_BaseModel):
    "Superclass for feature groups (which in turn define valid feature names)"
    prefix: FeaturePrefix

    def fmt_feature(self, rest: str) -> FeatureKey:
        return FeatureKey(f"{self.prefix}_{rest}")


class ColumnSpec(_BaseModel):
    feature_name: NonEmptyStr
    description: FeatureDesc


class VCFColumns(_BaseModel):
    "Columns for a vcf file"
    qual: ColumnSpec = ColumnSpec(
        feature_name="QUAL",
        description=FeatureDesc("The value of the QUAL column in the VCF file"),
    )
    filter: ColumnSpec = ColumnSpec(
        feature_name="FILTER",
        description=FeatureDesc("The value of the FILTER column in the VCF file"),
    )
    info: ColumnSpec = ColumnSpec(
        feature_name="INFO",
        description=FeatureDesc("The value of the INFO column in the VCF file"),
    )
    len: ColumnSpec = ColumnSpec(
        feature_name="indel_length",
        description=FeatureDesc(
            "The difference between the lengths of the ALT and REF columns in the VCF file."
        ),
    )


class VCFGroup(_FeatureGroup):
    "Feature and column names for VCF files"
    prefix: FeaturePrefix = "VCF"
    columns: VCFColumns = VCFColumns()

    def fmt_name(
        self: Self,
        f: Callable[[VCFColumns], ColumnSpec],
    ) -> tuple[FeatureKey, FeatureDesc]:
        c = f(self.columns)
        return (self.fmt_feature(c.feature_name), c.description)

    @property
    def str_features(self) -> dict[FeatureKey, FeatureDesc]:
        return dict(self.fmt_name(x) for x in [lambda x: x.filter, lambda x: x.info])

    def field_features(
        self,
        format_fields: set[str],
    ) -> dict[FeatureKey, FeatureDesc]:
        return {
            self.fmt_feature(f): FeatureDesc(
                f"The value of field '{f}' in the FORMAT/SAMPLE columns of the vcf"
            )
            for f in format_fields
        }

    def features(self, format_fields: set[str]) -> dict[FeatureKey, FeatureDesc]:
        return {
            **self.str_features,
            **dict(self.fmt_name(x) for x in [lambda x: x.qual, lambda x: x.len]),
            **self.field_features(format_fields),
        }


class MapSuffixes(_BaseModel):
    "Suffixes corresponding to low-map regions"
    low: NonEmptyStr = "difficult_100bp"
    high: NonEmptyStr = "difficult_250bp"


class MapGroup(_FeatureGroup):
    "Feature and column names for low-map dataframes"
    prefix: FeaturePrefix = "MAP"
    suffixes: MapSuffixes = MapSuffixes()

    @property
    def low(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.low)

    @property
    def high(self) -> FeatureKey:
        return self.fmt_feature(self.suffixes.high)

    @property
    def features(self) -> dict[FeatureKey, FeatureDesc]:
        return {
            self.low: FeatureDesc(
                (
                    "A low hard-to-map region is a 100bp sequence which "
                    "maps to at least one other region in the genome with "
                    "at most two substitutions and one gap. This is a binary "
                    "feature, so if a variant intersect with a low mappability "
                    "region, it gets a 1, otherwise 0."
                )
            ),
            self.high: FeatureDesc(
                (
                    "A high hard-to-map region is a 250bp sequence which "
                    "exactly maps to at least one other region in the genome. "
                    "This is a binary feature, so if a variant intersect with "
                    "a low mappability region, it gets a 1, otherwise 0."
                )
            ),
        }


class HomopolySuffixes(_BaseModel):
    "Suffixes corresponding to homopolymer regions"
    len: NonEmptyStr = "length"
    imp_frac: NonEmptyStr = "imperfect_frac"


class HomopolyGroup(_FeatureGroup):
    "Feature and column names for homopolymer dataframes"
    prefix: FeaturePrefix = "HOMOPOL"
    bases: set[Base] = set(b for b in Base)
    suffixes: HomopolySuffixes = HomopolySuffixes()

    def fmt_name(self, b: Base, f: Callable[[HomopolySuffixes], str]) -> FeatureKey:
        return self.fmt_feature(f"{b.value}_{f(self.suffixes)}")

    @property
    def features(self) -> dict[FeatureKey, FeatureDesc]:
        def fmt_len(b: Base) -> FeatureDesc:
            return FeatureDesc(
                (
                    f"The length of a {b.value} homopolymer."
                    "This includes all homopolymers of the same base as well "
                    f"as 'imperfect homopolymers' which have non-{b.value}s "
                    "separated by at least 4bp. Note that minimum length of "
                    "a homopolymer is 4bp."
                )
            )

        def fmt_imp_frac(b: Base) -> FeatureDesc:
            return FeatureDesc(
                f"The number of of non-{b.value}s over the length of the homopolymer."
            )

        fs = [(lambda x: x.len, fmt_len), (lambda x: x.imp_frac, fmt_imp_frac)]
        return {self.fmt_name(b, nf): df(b) for (nf, df), b in product(fs, self.bases)}


class RMSKGroup(_FeatureGroup):
    "Feature and column names for repeat masker dataframes"
    prefix: FeaturePrefix = "REPMASK"

    def fmt_name(
        self,
        f: RMSKFile,
        grp: str,
        fam: str | None,
    ) -> PandasColumn:
        assert grp in f.class_families, f"{grp} not a valid RMSK class"
        rest = maybe(grp, lambda f: f"{grp}_{fam}", fam)
        return PandasColumn(self.fmt_feature(f"{rest}_length"))

    def features(self, f: RMSKFile) -> dict[FeatureKey, FeatureDesc]:
        def fmt(cls: str, fam: str | None) -> tuple[FeatureKey, FeatureDesc]:
            k = FeatureKey(self.fmt_name(f, cls, fam))
            d = (
                FeatureDesc(f"The length of an intersecting region in class '{cls}'")
                if fam is not None
                else FeatureDesc(
                    f"The length of an intersecting region from family '{fam}' in class '{cls}'"
                )
            )
            return (k, d)

        return dict(
            [fmt(c, None) for c in f.class_families]
            + [fmt(c, f) for c, fs in f.class_families.items() for f in fs]
        )


DescribedFeature = tuple[FeatureKey, FeatureDesc]
DescribedColumn = tuple[PandasColumn, FeatureDesc]


class MergedFeatureGroup(_FeatureGroup, Generic[X]):
    "Superclass for feature group which supports bedtools merge"
    operations: set[BedMergeOp]
    description: NonEmptyStr = "segmental duplications"
    columns: X

    @property
    def count_feature(self) -> DescribedFeature:
        return (
            FeatureKey(self.fmt_feature("count")),
            FeatureDesc(
                f"The number of overlapping {self.description} regions used to make this feature."
            ),
        )

    def fmt_col(self, f: Callable[[X], ColumnSpec]) -> tuple[PandasColumn, FeatureDesc]:
        c = f(self.columns)
        return (PandasColumn(self.fmt_feature(c.feature_name)), c.description)

    def fmt_merged_feature(self, middle: str, op: BedMergeOp) -> FeatureKey:
        return FeatureKey(f"{middle}_{op.value}")

    # TODO generalize this so it can take any sequence rather than a hardcoded
    # dict?
    def ops_product(
        self,
        xs: dict[FeatureKey, FeatureDesc],
    ) -> dict[FeatureKey, FeatureDesc]:
        return {
            FeatureKey(f"{f}_{o.value}"): FeatureDesc(f"{d} ({o.value})")
            for (f, d), o in product(xs.items(), self.operations)
        }

    def merged_features(
        self,
        column_getters: list[Callable[[X], ColumnSpec]],
        other_columns: dict[FeatureKey, FeatureDesc] = {},
    ) -> dict[FeatureKey, FeatureDesc]:
        cs = dict(
            [
                (FeatureKey(self.fmt_feature(c.feature_name)), c.description)
                for c in [f(self.columns) for f in column_getters]
            ]
        )
        cnt = self.count_feature
        return {**self.ops_product({**cs, **other_columns}), cnt[0]: cnt[1]}


class SegDupsColumns(_BaseModel):
    "Columns corresponding to the superdups files"
    alignL: ColumnSpec = ColumnSpec(
        feature_name="size",
        description=FeatureDesc(
            (
                "The fraction (a float between 0 and 1) between this "
                "segmental duplication and another in a different part "
                "of the genome. Corresponds to 'alignL' in superdups database."
            )
        ),
    )
    fracMatchIndel: ColumnSpec = ColumnSpec(
        feature_name="identity",
        description=FeatureDesc(
            (
                "Spaces/positions in the alignment (positive integer). "
                "Corresponds to 'fracMatchIndel' in superdups database."
            )
        ),
    )


class SegDupsGroup(MergedFeatureGroup[SegDupsColumns]):
    "Feature and column names for segdups dataframes"
    prefix: FeaturePrefix = "SEGDUP"
    columns: SegDupsColumns = SegDupsColumns()
    operations: set[BedMergeOp] = {BedMergeOp.MIN, BedMergeOp.MAX, BedMergeOp.MEAN}
    description: NonEmptyStr = "segmental duplications"

    @property
    def features(self) -> dict[FeatureKey, FeatureDesc]:
        return self.merged_features([lambda x: x.alignL, lambda x: x.fracMatchIndel])


class TandemRepeatColumns(_BaseModel):
    "Columns corresponding to the simple_repeats files"
    period: ColumnSpec = ColumnSpec(
        feature_name="unit_size",
        description=FeatureDesc(
            "Length of the repeat unit. Corresponds to 'period' in TRF."
        ),
    )
    copyNum: ColumnSpec = ColumnSpec(
        feature_name="unit_copies",
        description=FeatureDesc(
            "Mean number of copies of the repeated unit. Corresponds to 'copyNum' in TRF."
        ),
    )
    perMatch: ColumnSpec = ColumnSpec(
        feature_name="identity",
        description=FeatureDesc(
            "Percentage match (integer between 0 and 100). Corresponds to 'perMatch' in TRF."
        ),
    )
    perIndel: ColumnSpec = ColumnSpec(
        feature_name="per_indel_mismatch",
        description=FeatureDesc(
            "Percentage INDEL (integer between 0 and 100). Corresponds to 'perIndel' in TRF."
        ),
    )
    score: ColumnSpec = ColumnSpec(
        feature_name="score",
        description=FeatureDesc(
            "Alignment score (integer with minimum of 50). Corresponds to 'score' in TRF."
        ),
    )


class TandemRepeatGroup(MergedFeatureGroup[TandemRepeatColumns]):
    "Feature and column names for segdups dataframes"
    prefix: FeaturePrefix = "TR"
    columns: TandemRepeatColumns = TandemRepeatColumns()
    operations: set[BedMergeOp] = {BedMergeOp.MIN, BedMergeOp.MAX, BedMergeOp.MEDIAN}

    def _base_name(self, bs: list[Base]) -> DescribedFeature:
        b = "".join(b.value for b in bs)
        f = self.fmt_feature(f"percent_{b}")
        d = f"The percentage of {b} in the repeat (integer from 0 to 100)."
        return (FeatureKey(f), FeatureDesc(d))

    @property
    def A(self) -> DescribedFeature:
        return self._base_name([Base.A])

    @property
    def T(self) -> DescribedFeature:
        return self._base_name([Base.T])

    @property
    def C(self) -> DescribedFeature:
        return self._base_name([Base.C])

    @property
    def G(self) -> DescribedFeature:
        return self._base_name([Base.G])

    @property
    def AT(self) -> DescribedFeature:
        return self._base_name([Base.A, Base.T])

    @property
    def GC(self) -> DescribedFeature:
        return self._base_name([Base.G, Base.C])

    @property
    def AG(self) -> DescribedFeature:
        return self._base_name([Base.A, Base.G])

    @property
    def CT(self) -> DescribedFeature:
        return self._base_name([Base.C, Base.T])

    @property
    def length(self) -> DescribedFeature:
        return (
            FeatureKey(self.fmt_feature("length")),
            FeatureDesc("The length of the this merged region of tandem repeats"),
        )

    @property
    def features(self) -> dict[FeatureKey, FeatureDesc]:
        return self.merged_features(
            [
                lambda x: x.period,
                lambda x: x.copyNum,
                lambda x: x.perMatch,
                lambda x: x.perIndel,
                lambda x: x.score,
            ],
            dict(
                [
                    self.A,
                    self.T,
                    self.G,
                    self.C,
                    self.AT,
                    self.AG,
                    self.CT,
                    self.GC,
                    self.length,
                ]
            ),
        )


class VariableGroup(_FeatureGroup):
    "Feature and column names for vcf manually-assigned variable dataframes"
    prefix: NonEmptyStr = "VAR"
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

    @property
    def all_keys(self) -> list[VarKey]:
        return list(self.continuous) + list(self.categorical)

    @property
    def features(self) -> dict[FeatureKey, FeatureDesc]:
        # return dict(self.fmt_feature(x) for x in self.all_keys())
        conts = [(k, v.description) for k, v in self.continuous.items()]
        cats = [(k, v.description) for k, v in self.categorical.items()]
        return {
            self.fmt_feature(f): FeatureDesc("<No description>") if d is None else d
            for f, d in conts + cats
        }

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


class FormatField(_BaseModel):
    """Means to parse a field in the FORMAT/SAMPLE columns of a VCF.

    Members:
    field_name: member of the FORMAT column to be parsed from SAMPLE
    field_missing: value to use if 'field_name' is not in FORMAT
    """

    field_name: str
    field_missing: str | None = None


FormatFields = dict[str, FormatField | str | None]


class UnlabeledVCFQuery(VCFFile):
    "A vcf to be used as the query for a model without labels."
    refset: RefsetKey
    variables: dict[VarKey, VarVal]
    max_ref: Annotated[int, Field(ge=0)] = 50
    max_alt: Annotated[int, Field(ge=0)] = 50
    format_fields: FormatFields = {
        "VAF": FormatField(field_name="VAF"),
        "DP": FormatField(field_name="DP"),
        "GT": FormatField(field_name="GT"),
        "GQ": FormatField(field_name="GQ"),
    }


class LabeledVCFQuery(UnlabeledVCFQuery):
    "A vcf to be used as the query for a model with labels (requires benchmark)."
    benchmark: BenchKey


VCFQuery = UnlabeledVCFQuery | LabeledVCFQuery


class FeatureDefs(_BaseModel):
    "Defines valid feature names to be specified in models"
    label: NonEmptyStr = "label"
    vcf: VCFGroup = VCFGroup()
    mappability: MapGroup = MapGroup()
    homopolymers: HomopolyGroup = HomopolyGroup()
    repeat_masker: RMSKGroup = RMSKGroup()
    segdups: SegDupsGroup = SegDupsGroup()
    tandem_repeats: TandemRepeatGroup = TandemRepeatGroup()
    variables: VariableGroup

    @property
    def label_name(self: Self) -> PandasColumn:
        return PandasColumn(self.label)

    @property
    def non_summary_cols(self) -> list[PandasColumn]:
        return [PandasColumn(c) for c in BED_COLS] + [
            PandasColumn(x) for x in self.vcf.str_features
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


def assert_member(x: X, xs: Sequence[X]) -> None:
    assert x in xs, f"'{x}' is not one of {xs}"


def flatten_features(fs: dict[FeatureKey, Feature]) -> list[FeatureKey]:
    return [k if v.alt_name is None else v.alt_name for k, v in fs.items()]


def assert_keypattern(x: str) -> None:
    assert re.fullmatch(KEY_CONSTR, x), f"key '{x}' does not match {KEY_CONSTR}"


LabeledQueryMap = dict[LabeledQueryKey, LabeledVCFQuery]
UnlabeledQueryMap = dict[UnlabeledQueryKey, UnlabeledVCFQuery]
RefMap = dict[RefKey, Reference]
RefsetMap = dict[RefsetKey, Refset]
ModelMap = dict[ModelKey, Model]


class StratoMod(_BaseModel):
    "Root config for the stratomod pipeline."

    paths: Paths
    tools: Tools
    references: RefMap
    feature_definitions: FeatureDefs
    reference_sets: RefsetMap
    labeled_queries: LabeledQueries
    unlabeled_queries: UnlabeledQueries
    models: ModelMap

    @classmethod
    def _lookup_vcfquery(
        cls,
        labeled: dict[LabeledQueryKey, LabeledVCFQuery],
        unlabeled: dict[UnlabeledQueryKey, UnlabeledVCFQuery],
        k: QueryKey,
    ) -> VCFQuery:
        try:
            return labeled[cast(LabeledQueryKey, k)]
        except KeyError:
            return unlabeled[cast(UnlabeledQueryKey, k)]

    @classmethod
    def _merge_features(
        cls,
        fs: FeatureDefs,
        lm: LabeledQueryMap,
        um: UnlabeledQueryMap,
        rm: RefMap,
        rsm: RefsetMap,
        ks: list[QueryKey],
    ) -> dict[FeatureKey, FeatureDesc]:
        def to_keys(
            f: Callable[[QueryKey], dict[FeatureKey, FeatureDesc]]
        ) -> dict[FeatureKey, FeatureDesc]:
            return dict(ChainMap(*[f(k) for k in ks]))

        def querykey_to_vcf(k: QueryKey) -> dict[FeatureKey, FeatureDesc]:
            x = set(cls._lookup_vcfquery(lm, um, k).format_fields)
            return fs.vcf.features(x)

        def querykey_to_rmsk(k: QueryKey) -> dict[FeatureKey, FeatureDesc]:
            rsk = cls._lookup_vcfquery(lm, um, k).refset
            rk = rsm[rsk].ref
            rmsk = rm[rk].feature_data.repeat_masker
            return fs.repeat_masker.features(rmsk)

        vcf_features = to_keys(querykey_to_vcf)
        rmsk_features = to_keys(querykey_to_rmsk)
        return {
            **vcf_features,
            **rmsk_features,
            **fs.mappability.features,
            **fs.homopolymers.features,
            **fs.segdups.features,
            **fs.tandem_repeats.features,
            **fs.variables.features,
        }

    @classmethod
    def _merge_feature_names(
        cls,
        fs: FeatureDefs,
        lm: LabeledQueryMap,
        um: UnlabeledQueryMap,
        rm: RefMap,
        rsm: RefsetMap,
        ks: list[QueryKey],
    ) -> set[FeatureKey]:
        return set(cls._merge_features(fs, lm, um, rm, rsm, ks))

    @validator("reference_sets", each_item=True)
    def refsets_have_valid_refkeys(
        cls: Type[Self],
        v: Refset,
        values: dict[str, Any],
    ) -> Refset:
        try:
            assert_member(v.ref, list(cast(RefMap, values["references"])))
            assert_keypattern(v.ref)
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
            assert_member(v.refset, list(cast(RefsetMap, values["reference_sets"])))
            assert_keypattern(v.refset)
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
            var_root = cast(FeatureDefs, values["feature_definitions"]).variables
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
            refsets = cast(RefsetMap, values["reference_sets"])
            refs = cast(RefMap, values["references"])
        except KeyError:
            pass
        else:
            ref_key = refsets[v.refset].ref
            ref_benchmarks = refs[ref_key].benchmarks
            assert_member(v.benchmark, list(ref_benchmarks))
            assert_keypattern(v.benchmark)
        return v

    @validator("unlabeled_queries")
    def input_keys_unique(
        cls: Type[Self],
        v: UnlabeledQueries,
        values: dict[str, Any],
    ) -> UnlabeledQueries:
        try:
            assert set(v).isdisjoint(
                set(cast(LabeledQueryMap, values["labeled_queries"]))
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
            query_keys = list(v.train) + [t.query_key for t in v.test.values()]
            features = cls._merge_feature_names(
                cast(FeatureDefs, values["feature_definitions"]),
                cast(LabeledQueryMap, values["labeled_queries"]),
                cast(UnlabeledQueryMap, values["unlabeled_queries"]),
                cast(RefMap, values["references"]),
                cast(RefsetMap, values["reference_sets"]),
                query_keys,
            )
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
            assert_subset(
                set(v.train),
                set(cast(LabeledQueryMap, values["labeled_queries"])),
            )
            for t in v.train:
                assert_keypattern(t)
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
            tests = [t.query_key for t in v.test.values()]
            ls = set(cast(LabeledQueryMap, values["labeled_queries"]))
            us = set(cast(UnlabeledQueryMap, values["unlabeled_queries"]))
            all_queries = ls | us
            assert_subset(set(tests), all_queries)
            for t in tests:
                assert_keypattern(t)
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
            var_root = cast(FeatureDefs, values["feature_definitions"]).variables
        except KeyError:
            pass
        else:
            varpairs = [
                (varname, varval)
                for t in v.test.values()
                for varname, varval in t.variables.items()
            ]
            assert set([p[0] for p in varpairs]) == set(
                var_root.all_keys
            ), "missing variables"
            for varname, varval in varpairs:
                var_root.validate_variable(varname, varval)
        return v

    @validator("models")
    def models_keys_are_valid(
        cls: Type[Self],
        v: dict[ModelKey, Model],
    ) -> dict[ModelKey, Model]:
        for k in v:
            assert_keypattern(k)
        return v

    # @validator("models", each_item=True)
    # def run_and_test_keys_are_valid(cls: Type[Self], v: Model) -> Model:
    #     for rk, r in v.runs.items():
    #         assert_keypattern(rk)
    #         for tk in r.test:
    #             assert_keypattern(tk)
    #     return v

    # various mapping functions

    def refsetkey_to_ref(self, key: RefsetKey) -> Reference:
        return self.references[self.refsetkey_to_refkey(key)]

    def refsetkey_to_refset(self, key: RefsetKey) -> Refset:
        return self.reference_sets[key]

    def _querykey_to_input(self, k: QueryKey) -> VCFQuery:
        return self._lookup_vcfquery(
            self.labeled_queries,
            self.unlabeled_queries,
            k,
        )

    def refsetkey_to_refkey(self, key: RefsetKey) -> RefKey:
        return self.refsetkey_to_refset(key).ref

    def refsetkey_to_chr_indices(self, key: RefsetKey) -> list[ChrIndex]:
        f = self.refsetkey_to_refset(key).chr_filter
        return sorted([x for x in ChrIndex] if len(f) == 0 else list(f))

    def refsetkey_to_chr_name_mapper(
        self,
        k: RefsetKey,
        p: ChrPrefix,
    ) -> dict[str, int]:
        cs = self.refsetkey_to_chr_indices(k)
        return {c.chr_name_full(p): c.value for c in cs}

    def refsetkey_to_chr_filter(
        self,
        get_prefix: Callable[[Reference], ChrPrefix],
        key: RefsetKey,
    ) -> ChrFilter:
        indices = self.refsetkey_to_chr_indices(key)
        prefix = get_prefix(self.refsetkey_to_ref(key))
        return ChrFilter(prefix, indices)

    def refsetkey_to_sdf_chr_filter(self, key: RefsetKey) -> list[str]:
        prefix, indices = self.refsetkey_to_chr_filter(lambda r: r.sdf.chr_prefix, key)
        return list(i.chr_name_full(prefix) for i in indices)

    def benchkey_to_vcf_chr_prefix(self, rkey: RefsetKey, bkey: BenchKey) -> str:
        return self.refsetkey_to_ref(rkey).benchmarks[bkey].vcf.chr_prefix

    def benchkey_to_bed_chr_prefix(self, rkey: RefsetKey, bkey: BenchKey) -> str:
        return self.refsetkey_to_ref(rkey).benchmarks[bkey].bed.params.chr_prefix

    def refkey_to_feature_data(self, key: RefKey) -> FeatureData:
        return self.references[key].feature_data

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
        return self.feature_definitions.variables.parse_vars(vs)

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
        tkey: TestKey,
    ) -> dict[FeatureKey, float]:
        test = self.models[mkey].test[tkey]
        qs = self.querykey_to_variables(test.query_key)
        return {**qs, **self.feature_definitions.variables.parse_vars(test.variables)}

    def modelkey_to_train_querykeys(self, mkey: ModelKey) -> list[LabeledQueryKey]:
        return [t for t in self.models[mkey].train]

    def modelkey_to_test_querykeys(self, mkey: ModelKey) -> list[QueryKey]:
        return [t.query_key for t in self.models[mkey].test.values()]

    def testkey_to_querykey(self, mkey: ModelKey, tkey: TestKey) -> QueryKey:
        return self.models[mkey].test[tkey].query_key

    def modelkey_to_feature_names(self, mk: ModelKey) -> set[FeatureKey]:
        """Return features for a given run.

        For the most part, the total set of defined features is constant across
        all runs. However, each query vcf can independently specify how to parse
        the pesky FORMAT column and a few other features that are
        reference-dependent. Due to this variability, we need to intersect all
        the field names together between all train and test keys.

        """
        train = self.modelkey_to_train_querykeys(mk)
        test = self.modelkey_to_test_querykeys(mk)
        query_key = train + test

        return self._merge_feature_names(
            self.feature_definitions,
            self.labeled_queries,
            self.unlabeled_queries,
            self.references,
            self.reference_sets,
            query_key,
        )

    # src getters

    def benchkey_to_vcf_src(self, rk: RefKey, bk: BenchKey) -> FileSrc:
        return self.references[rk].benchmarks[bk].vcf.src

    def benchkey_to_bed_src(self, rk: RefKey, bk: BenchKey) -> FileSrc:
        return self.references[rk].benchmarks[bk].bed.src

    def refkey_to_mappability_src(self, rk: RefKey, high: bool) -> FileSrc:
        return (
            self.references[rk].feature_data.mappability.high.src
            if high
            else self.references[rk].feature_data.mappability.low.src
        )

    def refkey_to_segdups_src(self, rk: RefKey) -> FileSrc:
        return self.references[rk].feature_data.segdups.src

    def refkey_to_simreps_src(self, rk: RefKey) -> FileSrc:
        return self.references[rk].feature_data.tandem_repeats.src

    def refkey_to_rmsk_src(self, rk: RefKey) -> FileSrc:
        return self.references[rk].feature_data.repeat_masker.src

    # paths

    @property
    def _log_dir(self) -> Path:
        return self.paths.results / "log"

    def _resource_or_log_dir(self, log: bool) -> Path:
        return self._log_dir / "resource" if log else self.paths.resources

    def _result_or_log_dir(self, log: bool) -> Path:
        return self._log_dir / "results" if log else self.paths.results

    def _labeled_dir(self, labeled: bool) -> Path:
        d, w = ("labeled", "l_query_key") if labeled else ("unlabeled", "ul_query_key")
        return Path(d) / get_wildcard(w)

    def _query_dir(self, labeled: bool, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "query"
            / get_wildcard("refset_key")
            / self._labeled_dir(labeled)
        )

    def query_src_dir(self, labeled: bool, log: bool) -> Path:
        return self._resource_or_log_dir(log) / "queries" / self._labeled_dir(labeled)

    def ref_src_dir(self, log: bool) -> Path:
        return self._resource_or_log_dir(log) / "reference" / get_wildcard("ref_key")

    def bench_src_dir(self, log: bool) -> Path:
        return self.ref_src_dir(log) / "bench" / get_wildcard("bench_key")

    def features_src_dir(self, log: bool) -> Path:
        return self.ref_src_dir(log) / "feature_data"

    def tool_src_dir(self, log: bool) -> Path:
        return self._resource_or_log_dir(log) / "tools"

    def tool_res_dir(self, log: bool) -> Path:
        return self._result_or_log_dir(log) / "tools"

    def bench_res_dir(self, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "bench"
            / get_wildcard("refset_key")
            / get_wildcard("bench_key")
        )

    def features_res_dir(self, log: bool) -> Path:
        return self._result_or_log_dir(log) / "features" / get_wildcard("refset_key")

    def query_prepare_res_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "prepare"

    def query_parsed_res_dir(self, labeled: bool, log: bool) -> Path:
        return self._query_dir(labeled, log) / "parsed"

    def vcfeval_res_dir(self, log: bool) -> Path:
        return self._query_dir(True, log) / "vcfeval"

    def refset_res_dir(self, log: bool) -> Path:
        return self._result_or_log_dir(log) / "references" / get_wildcard("refset_key")

    def annotated_res_dir(self, labeled: bool, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "annotated_variants"
            / self._labeled_dir(labeled)
        )

    def model_train_res_dir(self, log: bool) -> Path:
        return (
            self._result_or_log_dir(log)
            / "model"
            / wildcard_format("{}-{}", "model_key", "vartype_key")
        )

    def model_test_res_dir(self, labeled: bool, log: bool) -> Path:
        return (
            self.model_train_res_dir(log)
            / "test"
            / ("labeled" if labeled else "unlabeled")
            / get_wildcard("test_key")
        )

    # wildcard expansion

    def querykey_is_labeled(self, k: QueryKey) -> bool:
        return k in self.labeled_queries

    # @property
    # def _all_runs(self) -> list[ModelRun]:
    #     return [r for m in self.models.values() for r in m.runs.values()]

    @property
    def all_labeled_querykeys(self) -> set[LabeledQueryKey]:
        train_keys = [k for m in self.models.values() for k in m.train]
        test_keys = [
            LabeledQueryKey(k)
            for m in self.models.values()
            for t in m.test.values()
            if self.querykey_is_labeled(k := t.query_key)
        ]
        return set(train_keys + test_keys)

    @property
    def all_unlabeled_querykeys(self) -> set[UnlabeledQueryKey]:
        return set(
            [
                UnlabeledQueryKey(k)
                for m in self.models.values()
                for t in m.test.values()
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
            (ModelKeyCombo(model_key, vartype_key), model)
            for model_key, model in self.models.items()
            for vartype_key in model.vartypes
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
                vartype_key=map(lambda x: x.run_combo.vartype_key.value, key_set),
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
                vartype_key=map(lambda x: x.run_combo.vartype_key.value, key_set),
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
                vartype_key=map(lambda x: x.run_combo.vartype_key.value, key_set),
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
            vartype_key=map(lambda x: x.run_combo.vartype_key.value, train_set),
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
