from typing import Any, TextIO, Callable
import common.config as cfg
from common.io import with_gzip_maybe, setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


def dot_to_blank(s: str) -> str:
    return "" if s == "." else s


def none_to_blank(s: str | None) -> str:
    return "" if s is None else s


def write_row(
    fo: TextIO,
    chrom: str,
    start: str,
    end: str,
    qual: str,
    info: str,
    filt: str,
    indel_length: str,
    parse_fields: list[str],
    const_fields: list[str],
    label: str | None,
) -> None:
    const_cols = [chrom, start, end, qual, info, filt, indel_length]
    label_col = [] if label is None else [label]
    cols = [*const_cols, *parse_fields, *const_fields, *label_col]
    fo.write("\t".join(cols) + "\n")


def line_to_bed_row(
    fo: TextIO,
    ls: list[str],
    vcf: cfg.UnlabeledVCFQuery,
    vtk: cfg.VartypeKey,
    parse_fields: list[tuple[str, cfg.FormatField]],
    const_field_values: list[str],
    label: str | None,
) -> bool:
    # CHROM = 0
    # POS = 1
    # ID = 2
    # REF = 3
    # ALT = 4
    # QUAL = 5
    # FILTER = 6
    # INFO = 7
    # FORMAT = 8
    # SAMPLE = 9

    chrom = int(ls[0])
    start = int(ls[1]) - 1  # bed's are 0-indexed and vcf's are 1-indexed

    # remove cases where ref and alt are equal (which is what "." means)
    if ls[4] == "." or ls[3] == ls[4]:
        logger.info("Skipping equal variant at %s, %s", chrom, start)
        return False

    # remove multiallelics
    if "," in ls[4]:
        logger.info("Skipping multiallelic variant at %s, %s", chrom, start)
        return False

    # remove anything that doesn't pass out length filters
    ref_len = len(ls[3])
    alt_len = len(ls[4])

    if len(ls[3]) > vcf.max_ref or len(ls[4]) > vcf.max_alt:
        logger.info("Skipping oversized variant at %s, %s", chrom, start)
        return False

    # keep only the variant type we care about
    is_snv = ref_len == alt_len == 1

    if is_snv and vtk is cfg.VartypeKey.SNV:
        indel_length = 0
    elif not is_snv and ref_len != alt_len and vtk is cfg.VartypeKey.INDEL:
        indel_length = alt_len - ref_len
    else:
        return False

    # parse the format/sample columns if desired
    if len(parse_fields) > 0:
        fmt_col = ls[8].split(":")
        smpl_col = ls[9].split(":")
        # ASSUME any FORMAT/SAMPLE columns with different lengths are screwed
        # up in some way
        if len(fmt_col) != len(smpl_col):
            logger.error(
                "FORMAT/SAMPLE have different cardinality at %s, %s", chrom, start
            )
            return True
        d = dict(zip(fmt_col, smpl_col))
        parsed_field_values = [
            d[field.field_name]
            if field.field_name in d
            else none_to_blank(field.field_missing)
            for col, field in parse_fields
        ]
    else:
        parsed_field_values = []

    write_row(
        fo,
        str(chrom),
        str(start),
        str(start + ref_len),
        dot_to_blank(ls[5]),
        dot_to_blank(ls[6]),
        dot_to_blank(ls[7]),
        str(indel_length),
        parsed_field_values,
        list(const_field_values),
        label,
    )

    return False


def parse(smk: Any, sconf: cfg.StratoMod, fi: TextIO, fo: TextIO) -> None:
    defs = sconf.feature_definitions
    vcf = sconf.querykey_to_vcf(cfg.LabeledQueryKey(smk.params.query_key))
    vtk = cfg.VartypeKey(smk.wildcards.vartype_key)
    found_error = False

    try:
        label = str(smk.wildcards.label)
    except AttributeError:
        label = None

    fields = [(str(defs.vcf.fmt_feature(k)), v) for k, v in vcf.format_fields.items()]
    parse_fields = [(k, v) for k, v in fields if isinstance(v, cfg.FormatField)]
    const_fields = [
        (k, none_to_blank(v)) for k, v in fields if not isinstance(v, cfg.FormatField)
    ]
    # unzip only works on non-empty lists :(
    const_field_names = [f[0] for f in const_fields]
    const_field_values = [f[1] for f in const_fields]

    # write header
    def fmt(f: Callable[[cfg.VCFColumns], cfg.ColumnSpec]) -> str:
        return defs.vcf.fmt_name(f)[0]

    write_row(
        fo,
        cfg.BED_CHROM,
        cfg.BED_START,
        cfg.BED_END,
        fmt(lambda x: x.qual),
        fmt(lambda x: x.info),
        fmt(lambda x: x.filter),
        fmt(lambda x: x.len),
        [f[0] for f in parse_fields],
        list(const_field_names),
        None if label is None else defs.label_name,
    )

    for ln in fi:
        if ln.startswith("#"):
            continue

        err = line_to_bed_row(
            fo,
            ln.rstrip().split("\t"),
            vcf,
            vtk,
            parse_fields,
            list(const_field_values),
            label,
        )

        found_error = err or found_error

    if found_error is True:
        exit(1)


def main(smk: Any, config: cfg.StratoMod) -> None:
    with_gzip_maybe(
        lambda i, o: parse(smk, config, i, o),
        str(smk.input[0]),
        str(smk.output[0]),
    )


main(snakemake, snakemake.config)  # type: ignore