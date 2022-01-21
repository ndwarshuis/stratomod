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
    "CHROM\tPOS\tPOS+length(REF)\tFILTER\tGT\tGQ\tDP\tVAF\tindel_length\tlabel\n"
)
f_out.flush()

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
    filter = split_line[6]
    sample = split_line[9]
    sample_split = sample.split(":")
    GT = sample_split[0]
    GQ = sample_split[1]
    DP = sample_split[2]
    VAF = sample_split[4]
    pos_plus_length_ref = int(pos) + len(alt)
    to_write_out = (
        chrom
        + "\t"
        + str(pos)
        + "\t"
        + str(pos_plus_length_ref)
        + "\t"
        + filter
        + "\t"
        + GT
        + "\t"
        + GQ
        + "\t"
        + DP
        + "\t"
        + VAF
        + "\t"
        + str(indel_length)
        + "\t"
        + args.label
        + "\n"
    )
    f_out.write(to_write_out)
    f_out.flush()


f.close()
f_out.close()
