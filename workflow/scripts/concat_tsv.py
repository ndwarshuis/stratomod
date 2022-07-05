import pandas as pd
from common.config import fmt_vcf_feature
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically


vaf_col = fmt_vcf_feature(snakemake.config, "vaf")
ad_col = fmt_vcf_feature(snakemake.config, "ad")
dp_col = fmt_vcf_feature(snakemake.config, "dp")


# ASSUME when using AD to calculate DP and VAF that we have already filtered
# out multi-allelic entries, and therefore these should all be scalers and not
# (nested) vectors
def ad_split(row):
    return [*map(float, row[ad_col].split(","))]


def ad_to_vaf(row):
    s = ad_split(row)
    return s[1] / (s[0] + s[1])


def ad_to_dp(row):
    s = ad_split(row)
    return s[0] + s[1]


def fill_vaf_dp(df):
    # this may not all be necessary, but I want to make sure I can trace back
    # VAF/DP in case the calculated and non-calculated versions don't match
    # for whatever reason
    vaf_nil = df[vaf_col].isnull()
    dp_nil = df[dp_col].isnull()
    ad_nil = df[ad_col].isnull()
    vaf_mask = vaf_nil & ~ad_nil
    dp_mask = dp_nil & ~ad_nil
    df.loc[vaf_mask, vaf_col] = df[vaf_mask].apply(ad_to_vaf, axis=1)
    df.loc[vaf_mask, dp_col] = df[vaf_mask].apply(ad_to_dp, axis=1)
    df["_VAF_filled"] = vaf_mask
    df["_DP_filled"] = dp_mask
    return df


df = pd.concat([read_tsv(i, header=0) for i in snakemake.input])
write_tsv(snakemake.output[0], sort_bed_numerically(fill_vaf_dp(df)))
