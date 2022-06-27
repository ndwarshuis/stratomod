import pandas as pd
from common.config import fmt_vcf_feature
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically


vaf_col = fmt_vcf_feature(snakemake.config, "vaf")
ad_col = fmt_vcf_feature(snakemake.config, "ad_sum")
dp_col = fmt_vcf_feature(snakemake.config, "dp")


def fill_vaf_dp(df):
    # this may not all be necessary, but I want to make sure I can trace back
    # VAF/DP in case the calculated and non-calculated versions don't match
    # for whatever reason
    vaf_nil = df[vaf_col].isnull()
    dp_nil = df[dp_col].isnull()
    ad_nil = df[ad_col].isnull()
    vaf_mask = vaf_nil & ~dp_nil & ~ad_nil
    dp_mask = ~vaf_nil & dp_nil & ~ad_nil
    df.loc[vaf_mask, vaf_col] = df.loc[vaf_mask, ad_col] / df.loc[vaf_mask, dp_col]
    df.loc[dp_mask, dp_col] = df.loc[dp_mask, ad_col] / df.loc[dp_mask, vaf_col]
    df["_VAF_filled"] = vaf_mask
    df["_AD_filled"] = dp_mask
    return df


df = pd.concat([read_tsv(i, header=0) for i in snakemake.input])
write_tsv(snakemake.output[0], sort_bed_numerically(fill_vaf_dp(df)))
