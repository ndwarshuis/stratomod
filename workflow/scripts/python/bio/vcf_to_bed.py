from typing import Any, TextIO
import common.config as cfg
from common.io import with_gzip_maybe, setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


def is_real(s: str) -> bool:
    return s.removeprefix("-").replace(".", "", 1).isdigit()


def dot_to_blank(s: str) -> str:
    return "" if s == "." else s


def none_to_blank(s: str | None) -> str:
    return "" if s is None else s


def write_row(
    fo: TextIO,
    chrom: str,
    start: str,
    end: str,
    vid: str,
    ref: str,
    alt: str,
    qual: str,
    filt: str,
    info: str,
    indel_length: str,
    parse_fields: list[str],
    const_fields: list[str],
    label: str | None,
) -> None:
    const_cols = [chrom, start, end, vid, ref, alt, qual, info, filt, indel_length]
    label_col = [] if label is None else [label]
    cols = [*const_cols, *parse_fields, *const_fields, *label_col]
    fo.write("\t".join(cols) + "\n")


def lookup_field(f: cfg.FormatField, d: dict[str, str]) -> str:
    try:
        v = d[f.name]
        if len(f.mapper) == 0:
            return v if is_real(v) else ""
        try:
            return str(f.mapper[v])
        except KeyError:
            return ""
    except KeyError:
        return none_to_blank(f.missing)


def line_to_bed_row(
    fo: TextIO,
    ls: list[str],
    vcf: cfg.UnlabeledVCFQuery,
    vtk: cfg.VartypeKey,
    parse_fields: list[cfg.FormatField],
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
    ref = ls[3]
    alt = ls[4]

    # remove cases where ref and alt are equal (which is what "." means)
    if alt == "." or ref == alt:
        logger.info("Skipping equal variant at %s, %s", chrom, start)
        return False

    # remove multiallelics
    if "," in alt:
        logger.info("Skipping multiallelic variant at %s, %s", chrom, start)
        return False

    # remove anything that doesn't pass out length filters
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len > vcf.max_ref or alt_len > vcf.max_alt:
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
            logger.error("FORMAT/SAMPLE have different lengths at %s, %s", chrom, start)
            return True
        d = dict(zip(fmt_col, smpl_col))
        parsed_field_values = [lookup_field(f, d) for f in parse_fields]
    else:
        parsed_field_values = []

    write_row(
        fo,
        str(chrom),
        str(start),
        str(start + ref_len),
        dot_to_blank(ls[2]),
        dot_to_blank(ref),
        dot_to_blank(alt),
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
        lbl = cfg.VCFLabel(smk.wildcards.label)
        if lbl is cfg.VCFLabel.TPBL:
            lbl = cfg.VCFLabel.TP
        label = lbl.value
    except AttributeError:
        label = None

    fields = [(str(defs.vcf.fmt_feature(k)), v) for k, v in vcf.format_fields.items()]
    parse_fields = [(k, v) for k, v in fields if isinstance(v, cfg.FormatField)]
    const_fields = [
        (k, none_to_blank(v)) for k, v in fields if not isinstance(v, cfg.FormatField)
    ]

    # write header
    write_row(
        fo,
        cfg.BED_CHROM,
        cfg.BED_START,
        cfg.BED_END,
        defs.vcf.id,
        defs.vcf.ref,
        defs.vcf.alt,
        defs.vcf.qual[0],
        defs.vcf.filter,
        defs.vcf.info,
        defs.vcf.indel_length[0],
        [f[0] for f in parse_fields],
        [f[0] for f in const_fields],
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
            [f[1] for f in parse_fields],
            [f[1] for f in const_fields],
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
