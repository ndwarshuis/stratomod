from functools import partial
from collections import namedtuple
import numpy as np
import pandas as pd
from more_itertools import partition
from common.tsv import write_tsv
from common.bed import standardize_chr_series
from common.cli import setup_logging
from common.config import fmt_vcf_feature, assert_empty, lookup_train_test_input
from common.prepare import compose

logger = setup_logging(snakemake.log[0])

input_col = namedtuple("InputCol", ["dtype", "na_value"])


def read_vcf(input_cols, path):
    ecs = [*enumerate(input_cols.items())]
    dtypes = dict((i, c.dtype) for i, (_, c) in ecs)
    na_values = dict((i, c.na_value) for i, (_, c) in ecs if c.na_value is not None)
    # Add columns manually because pandas can't distinguish between lines
    # starting with '##' or '#'
    df = pd.read_table(
        path,
        comment="#",
        header=None,
        na_values=na_values,
        dtype=dtypes,
        names=[*input_cols],
    )
    # NOTE: if FORMAT/SAMPLE don't exist these will just be NaNs
    return df


def standardize_chrs(chr_col, df):
    df[chr_col] = standardize_chr_series(df[chr_col])
    return df.dropna(subset=[chr_col]).astype({chr_col: int})


def fix_dot_alts(df):
    # TODO not sure if these should be filled or outright removed
    df["ALT"] = df["REF"].where(df["ALT"].isna(), df["ALT"])
    return df


def remove_multialleles(df):
    # ASSUME there are no NaNs in alt at this point
    return df[~df["ALT"].str.contains(",")]


def lookup_maybe(k, d):
    return d[k] if k in d else np.nan


def assign_format_sample_fields(sconf, chrom, start, fields, df):
    fill_fields_gen, nan_fields_gen = partition(
        lambda x: x[1] is None,
        ((fmt_vcf_feature(sconf, k), v) for k, v in fields.items()),
    )
    for dest_col, _ in nan_fields_gen:
        df[dest_col] = np.nan

    fill_fields = [*fill_fields_gen]

    # Any fields that we don't want to parse will just get filled with NaNs. If
    # all fields are like this, our job is super easy and we don't need to
    # bother with zipping the FORMAT/SAMPLE columns
    if len(fill_fields) == 0:
        return df

    # split format and sample fields into lists
    def split_col(n):
        split_ser = df[n].str.split(":")
        len_ser = split_ser.str.len()
        return split_ser, len_ser

    format_split, format_len = split_col("FORMAT")
    sample_split, sample_len = split_col("SAMPLE")

    # warn user if any FORMAT/SAMPLE pairs are different length (which
    # shouldn't happen, but vcf files are weird)
    present = ~format_len.isna() & ~sample_len.isna()
    eqlen = format_len == sample_len

    assert_empty(
        [
            f"{r[chrom]}@{r[start]}"
            for r in df[~eqlen & present].to_dict(orient="records")
        ],
        "Lines where FORMAT/SAMPLE have different cardinality",
    )

    # zip FORMAT/SAMPLE into dict series and assign
    #
    # NOTE: this is a faily slow operation because we need to (or at least
    # should) not assume that FORMAT/SAMPLE pair is exactly the same wrt
    # order, cardinality, etc. Some fields that are in one line might be
    # absent from others. VCF files are weird...
    mask = eqlen & present
    splits = format_split[mask].combine(
        sample_split[mask],
        lambda a, b: dict(zip(a, b)),
    )
    return pd.concat(
        [df]
        + [
            splits.apply(partial(lookup_maybe, field)).rename(dest_col)
            for dest_col, field in fill_fields
        ],
        axis=1,
    )


def filter_mask(ref_len, alt_len, filter_key):
    assert filter_key in ["INDEL", "SNP"], f"Invalid filter key: {filter_key}"
    snps = (ref_len == 1) & (alt_len == 1)
    return ~snps & ~(ref_len == alt_len) if filter_key == "INDEL" else snps


def add_length_and_filter(indel_len_col, start_col, end_col, filter_key, df):
    # ASSUME there are no NaNs in alt at this point
    multi_alts = df["ALT"].str.contains(",")
    alt_len = df["ALT"].str.len()
    ref_len = df["REF"].str.len()
    df[indel_len_col] = alt_len - alt_len
    df[end_col] = df[start_col] + ref_len
    mask = filter_mask(ref_len, alt_len, filter_key) & ~multi_alts
    return df[mask].copy()


def select_nonlabel_columns(sconf, non_field_cols, fields, df):
    split_keep = [fmt_vcf_feature(sconf, f) for f in fields]
    return df[non_field_cols + split_keep]


def add_label_maybe(wildcards, label_col, df):
    if hasattr(wildcards, "label"):
        df[label_col] = wildcards.label
    return df


def main():
    wildcards = snakemake.wildcards
    sconf = snakemake.config
    fconf = sconf["features"]
    iconf = lookup_train_test_input(sconf, wildcards.input_key)
    idx = fconf["index"]

    chrom = idx["chr"]
    pos = idx["start"]
    qual = fmt_vcf_feature(sconf, "qual")
    filt = fmt_vcf_feature(sconf, "filter")
    end = idx["end"]
    indel_length = fmt_vcf_feature(sconf, "len")

    input_cols = {
        chrom: input_col(str, None),
        pos: input_col(int, None),
        "ID": input_col(str, "."),
        "REF": input_col(str, None),
        "ALT": input_col(str, "."),
        qual: input_col(float, "."),
        filt: input_col(str, "."),
        "INFO": input_col(str, "."),
        "FORMAT": input_col(str, "."),
        "SAMPLE": input_col(str, "."),
    }

    non_field_cols = [chrom, pos, end, indel_length, qual, filt]

    fields = iconf["format_fields"]

    return compose(
        partial(write_tsv, snakemake.output[0]),
        partial(add_label_maybe, wildcards, fconf["label"]),
        partial(select_nonlabel_columns, sconf, non_field_cols, fields),
        partial(assign_format_sample_fields, sconf, chrom, pos, fields),
        partial(standardize_chrs, chrom),
        partial(
            add_length_and_filter,
            fmt_vcf_feature(sconf, "len"),
            pos,
            end,
            wildcards.filter_key,
        ),
        remove_multialleles,
        fix_dot_alts,
        partial(read_vcf, input_cols),
    )(snakemake.input[0])


main()


# header = make_header(label_val)

# f = open(snakemake.input[0], "r")
# f_out = open(snakemake.output[0], "w+")
# lines = f.readlines()

# f_out.write("{}\n".format("\t".join(header)))
# f_out.flush()


# # def const_na(_):
# #     return NAN


# parse = lookup_train_test_input(sconf, wcs.input_key)["parse"]

# vaf_parser = (
#     partial(lookup_maybe, "VAF") if parse is not None and parse["vaf"] else const_na
# )
# dp_parser = (
#     partial(lookup_maybe, "DP") if parse is not None and parse["dp"] else const_na
# )


# # TODO different VCFs have different fields, we want to have DP and VAF almost
# # always, can (almost always) just use the AD field to get the VAF
# for line in lines:
#     if line.startswith("#"):
#         continue
#     split_line = line.split("\t")
#     chrom = split_line[0]
#     pos = split_line[1]
#     start = 0
#     end = 0
#     ref = split_line[3]
#     alt = split_line[4]
#     qual = split_line[5]
#     # CHROM, POS, POS+length(REF), FILTER, GT, GQ, DP, and VAF (only DP and VAF will probably be used inputs to the EBM - the rest are for our info)
#     ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002
#     # chr1    631859  .       CG      C       46.8    PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:41:34:1,33:0.970588:46,41,0
#     if "," in alt:
#         continue
#     ref_length = len(ref)
#     alt_length = len(alt)
#     # if we want INDELs skip everything that has REF/ALT of one BP or the same
#     # number of BPs
#     if wcs.filter_key == "INDEL" and (
#         (ref_length == alt_length == 1) or ref_length == alt_length
#     ):
#         continue
#     # if we want SNPs, skip everything that isn't REF/ALT with one BP
#     if wcs.filter_key == "SNP" and not (ref_length == alt_length == 1):
#         continue
#     indel_length = alt_length - ref_length
#     filt = split_line[6]
#     if parse is not None:
#         fmt = split_line[8].split(":")
#         # rstrip the newline off at the end
#         sample = split_line[9].rstrip().split(":")
#         if len(fmt) != len(sample):
#             logger.warn(
#                 "FORMAT/SAMPLE have different cardinality: %s %d %d",
#                 chrom,
#                 start,
#                 end,
#             )
#             continue

#         named_sample = dict(zip(fmt, sample))
#     else:
#         named_sample = {}

#     pos_plus_length_ref = int(pos) + len(alt)
#     to_write_out = "\t".join(
#         [
#             chrom,
#             str(pos),
#             str(pos_plus_length_ref),
#             qual,
#             filt,
#             lookup_maybe("GT", named_sample),
#             lookup_maybe("GQ", named_sample),
#             dp_parser(named_sample),
#             vaf_parser(named_sample),
#             str(indel_length),
#             *label_val,
#         ]
#     )
#     f_out.write(f"{to_write_out}\n")
#     f_out.flush()


# f.close()
# f_out.close()
