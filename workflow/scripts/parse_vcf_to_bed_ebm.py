from functools import partial
from common.cli import setup_logging
from common.config import fmt_vcf_feature

logger = setup_logging(snakemake.log[0])

wcs = snakemake.wildcards
has_label = hasattr(wcs, "label")
label_val = [wcs.label] if has_label else []
vartype = wcs.filter_key

sconf = snakemake.config
fconf = sconf["features"]
vconf = fconf["vcf"]
idx = fconf["index"]


def make_header(has_label):
    base_cols = [
        idx["chr"],
        idx["start"],
        idx["end"],
        *map(
            lambda f: fmt_vcf_feature(sconf, f),
            ["qual", "filter", "gt", "gq", "dp", "vaf", "len"],
        ),
    ]
    return [*base_cols, fconf["label"]] if has_label else base_cols


header = make_header(label_val)

f = open(snakemake.input[0], "r")
f_out = open(snakemake.output[0], "w+")
lines = f.readlines()

f_out.write("{}\n".format("\t".join(header)))
f_out.flush()

NAN = "NaN"


def lookup_maybe(k, d):
    return d[k] if k in d else NAN


def const_na(_):
    return NAN


parse = sconf["inputs"][wildcards.input_key]["parse"]

vaf_parser = partial(lookup_maybe, "VAF") if parse["vaf"] else const_na
dp_parser = partial(lookup_maybe, "DP") if parse["dp"] else const_na


# TODO different VCFs have different fields, we want to have DP and VAF almost
# always, can (almost always) just use the AD field to get the VAF
for line in lines:
    if line.startswith("#"):
        continue
    split_line = line.split("\t")
    chrom = split_line[0]
    pos = split_line[1]
    start = 0
    end = 0
    ref = split_line[3]
    alt = split_line[4]
    qual = split_line[5]
    # CHROM, POS, POS+length(REF), FILTER, GT, GQ, DP, and VAF (only DP and VAF will probably be used inputs to the EBM - the rest are for our info)
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002
    # chr1    631859  .       CG      C       46.8    PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:41:34:1,33:0.970588:46,41,0
    if "," in alt:
        continue
    ref_length = len(ref)
    alt_length = len(alt)
    # if we want INDELs skip everything that has REF/ALT of one BP or the same
    # number of BPs
    if wildcards.filter_key == "INDEL" and (
        (ref_length == alt_length == 1) or ref_length == alt_length
    ):
        continue
    # if we want SNPs, skip everything that isn't REF/ALT with one BP
    if wildcards.filter_key == "SNP" and not (ref_length == alt_length == 1):
        continue
    indel_length = alt_length - ref_length
    filt = split_line[6]
    fmt = split_line[8].split(":")
    # rstrip the newline off at the end
    sample = split_line[9].rstrip().split(":")
    if len(fmt) != len(sample):
        logger.warn(
            "FORMAT/SAMPLE have different cardinality: %s %d %d",
            chrom,
            start,
            end,
        )
        continue

    named_sample = dict(zip(fmt, sample))
    pos_plus_length_ref = int(pos) + len(alt)
    to_write_out = "\t".join(
        [
            chrom,
            str(pos),
            str(pos_plus_length_ref),
            qual,
            filt,
            lookup_maybe("GT", named_sample),
            lookup_maybe("GQ", named_sample),
            dp_parser(named_sample),
            vaf_parser(named_sample),
            str(indel_length),
            *label_val,
        ]
    )
    f_out.write(f"{to_write_out}\n")
    f_out.flush()


f.close()
f_out.close()
