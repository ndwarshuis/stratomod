from typing import List, Dict, Optional, Set, Union, NewType
from enum import Enum
from pydantic import BaseModel

RefKey = NewType("RefKey", str)
RefsetKey = NewType("RefsetKey", str)
TestKey = NewType("TestKey", str)
InputKey = NewType("InputKey", str)
BenchKey = NewType("BenchKey", str)
ModelKey = NewType("ModelKey", str)
RunKey = NewType("RunKey", str)

VarKey = NewType("VarKey", str)


class Bases(Enum):
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


class EBMRun(BaseModel):
    inputs: Dict[RunKey, ModelInput]
    filter: Set[FilterKey]
    ebm_settings: EBMSettings
    error_labels: Set[Label]
    filtered_are_candidates: bool
    interactions: List[Union[str, List[str]]]
    features: Dict[str, Feature]


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


class MapColumns(BaseModel):
    low: str
    high: str


class MapMeta(BaseModel):
    prefix: str
    columns: MapColumns


class HomopolySuffixes(BaseModel):
    len: str
    imp_frac: str


class HomopolyMeta(BaseModel):
    prefix: str
    bases: Bases
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
