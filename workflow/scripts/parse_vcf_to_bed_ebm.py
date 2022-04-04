import argparse

parser = argparse.ArgumentParser(description="Parse VCF to BED")
parser.add_argument("--input", metavar="I", type=str, nargs="+", help="input vcf file")
parser.add_argument(
    "--output", metavar="I", type=str, nargs="+", help="output bed file"
)
parser.add_argument("--type", metavar="t", type=str, help="'INDEL' or 'SNP'")
parser.add_argument("--label", metavar="L", type=str, help="the label to assign")
args = parser.parse_args()

f = open(args.input[0], "r")
f_out = open(args.output[0], "w+")

lines = f.readlines()

f_out.write(
    "{}\n".format(
        "\t".join(
            [
                "CHROM",
                "POS",
                "POS+length(REF)",
                "QUAL",
                "FILTER",
                "GT",
                "GQ",
                "DP",
                "VAF",
                "indel_length",
                "label",
            ]
        )
    )
)
f_out.flush()


def lookup_maybe(d, k):
    return d[k] if k in d else "NaN"


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
    if args.type == "INDEL" and (
        (ref_length == alt_length == 1) or ref_length == alt_length
    ):
        continue
    # if we want SNPs, skip everything that isn't REF/ALT with one BP
    if args.type == "SNP" and not (ref_length == alt_length == 1):
        continue
    indel_length = alt_length - ref_length
    filt = split_line[6]
    fmt = split_line[8].split(":")
    # rstrip the newline off at the end
    sample = split_line[9].rstrip().split(":")
    if len(fmt) != len(sample):
        print(
            "WARNING: FORMAT/SAMPLE have different cardinality: {} {} {}".format(
                chrom, start, end
            )
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
            lookup_maybe(named_sample, "GT"),
            lookup_maybe(named_sample, "GQ"),
            lookup_maybe(named_sample, "DP"),
            lookup_maybe(named_sample, "VAF"),
            str(indel_length),
            args.label,
        ]
    )
    f_out.write(f"{to_write_out}\n")
    f_out.flush()


f.close()
f_out.close()
