# commands for initial EBM VCF parsing
# Done using bedtools v2.27.1

# conda install -c bioconda rtg-tools

# Intermediate files available at https://docs.google.com/document/d/1PNqFuqymoe9kduCbPlXJ_DEvtukEv2SpioohOgZEDLQ/edit 
# /aigenomics/deepvariant_output/HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22.vcf.gz

```
sed -e '/.RefCall./ s/\.\/\./0\/1/g' HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22.vcf |  sed -e '/.RefCall./ s/0\/0/0\/1/g' > HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22_ready_for_vcfeval.vcf
```
# GRCh38.sdf from https://s3.amazonaws.com/rtg-datasets/references/GRCh38.sdf.zip
# HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz and HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/

```
rtg vcfeval -b HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22_ready_for_vcfeval.vcf.gz -o truth_HG002_v4.2.1_query_HG002_DV_Illumina -t GRCh38.sdf --ref-overlap --all-records

python parse_vcf_to_bed_ebm_tp.py --input tp.vcf --output HG002_DV_Illumina_tp_ebm.bed

python parse_vcf_to_bed_ebm_fp.py --input fp.vcf --output HG002_DV_Illumina_fp_ebm.bed

cp HG002_DV_Illumina_tp_ebm.bed HG002_DV_Illumina_df_ebm.bed

cat HG002_DV_Illumina_fp_ebm.bed >> HG002_DV_Illumina_df_ebm.bed

cut -f7-9 HG002_DV_Illumina_df_ebm.bed > HG002_DV_Illumina_df_ebm_DP_VAF_label.bed
```

# Simple training test - run in interactive session in conda environment after running `pip install interpret`
```
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

df = pd.read_csv("HG002_DV_Illumina_df_ebm_DP_VAF_label.bed", sep="\t")
df.columns = ["DP", "VAF", "label"]

df = df.sample(frac=0.05)
train_cols = df.columns[0:2]
label = df.columns[-1]
X = df[train_cols]
y = df[label]

seed = 1

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=seed)

ebm = ExplainableBoostingClassifier(random_state=seed)
ebm.fit(X_train, y_train)

ebm_global = ebm.explain_global()
show(ebm_global)

ebm_local = ebm.explain_local(X_test[:5], y_test[:5])
show(ebm_local)

ebm.predict_and_contrib(X_test, output='probabilities')
```



# Annotate with genomicSuperDups.txt downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

```
cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 4 -o min,max,count,mean > genomicSuperDups_alignL_stats.bed

cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 5 -o min,max,count,mean > genomicSuperDups_fracMatchIndel_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,13,14,15,16,21,22,23,24 > HG002_DV_Illumina_df_ebm_annotated_with_segdups.bed
```