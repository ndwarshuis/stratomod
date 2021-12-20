# EBM VCF parsing and generating dataframe for training

## dependencies

 - bedtools => v2.27.1
 - conda install -c bioconda rtg-tools

# Intermediate files available at https://docs.google.com/document/d/1PNqFuqymoe9kduCbPlXJ_DEvtukEv2SpioohOgZEDLQ/edit 
# /aigenomics/deepvariant_output/HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22.vcf.gz

```
sed -e '/.RefCall./ s/\.\/\./0\/1/g' HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22.vcf |  sed -e '/.RefCall./ s/0\/0/0\/1/g' > HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22_ready_for_vcfeval.vcf
```
# GRCh38.sdf from https://s3.amazonaws.com/rtg-datasets/references/GRCh38.sdf.zip
# HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz and HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/

# install https://anaconda.org/bioconda/rtg-tools

```
rtg vcfeval -b HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -e HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c HG002.hiseqx.pcr-free.40x.dedup.grch38_chr1_22_ready_for_vcfeval.vcf.gz -o truth_HG002_v4.2.1_query_HG002_DV_Illumina -t GRCh38.sdf --ref-overlap --all-records

python parse_vcf_to_bed_ebm_tp.py --input tp.vcf --output HG002_DV_Illumina_tp_ebm_snps.bed

python parse_vcf_to_bed_ebm_fp.py --input fp.vcf --output HG002_DV_Illumina_fp_ebm_snps.bed

cp HG002_DV_Illumina_tp_ebm_snps.bed HG002_DV_Illumina_df_ebm_snps.bed

cat HG002_DV_Illumina_fp_ebm_snps.bed >> HG002_DV_Illumina_df_ebm_snps.bed

cut -f7-9 HG002_DV_Illumina_df_ebm_snps.bed > HG002_DV_Illumina_df_ebm_snps_DP_VAF_label.bed
```



# Annotate with genomicSuperDups.txt downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

```
cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 4 -o min,max,count,mean > genomicSuperDups_alignL_stats.bed

cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 5 -o min,max,count,mean > genomicSuperDups_fracMatchIndel_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,13,14,15,16,21,22,23,24 > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups.bed
```

## EBM for SNPs
``` bash
bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps.bed -b GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed -b GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed

cat HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed | cut -f7,8,9,13,18 > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers.bed

HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers.bed

cut -f1,2 HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_DP_VAF.bed
cut -f3 HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_label.bed
cut -f4,5 HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_CG_AT_lens.bed
paste HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_DP_VAF.bed HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_CG_AT_lens.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_DP_VAF_CG_AT_lens_temp.bed
paste HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_DP_VAF_CG_AT_lens_temp.bed HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_label.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_homopolymers_DP_VAF_CG_AT_lens_label.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_snps_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,13,18,24,25,32 > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers.bed

cut -f1,2 HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF.bed
cut -f3 HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_label.bed
cut -f4,5,6,7,8 HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_CG_AT_lens.bed
paste HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF.bed HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_CG_AT_lens.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF_CG_AT_lens_temp.bed
paste HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF_CG_AT_lens_temp.bed HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_label.bed > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF_CG_AT_lens_label.bed
```


``` python
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

df = pd.read_csv("HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF_CG_AT_lens_alignLMax_alignLCount_fracMatchIndelMax_label.bed", sep="\t")
df.columns = ["DP", "VAF", "CGhomopolgt3_len", "AThomopolgt3_len", "alignL_max", "alignL_count", "fracMatchIndel_max", "label"]

# adapted from #3 at https://towardsdatascience.com/converting-data-to-a-numeric-type-in-pandas-db9415caab0b
df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype('int')
df['VAF'] = pd.to_numeric(df['VAF'], errors='coerce').fillna(0.0).astype('float')
df['CGhomopolgt3_len'] = pd.to_numeric(df['CGhomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['AThomopolgt3_len'] = pd.to_numeric(df['AThomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['alignL_max'] = pd.to_numeric(df['alignL_max'], errors='coerce').fillna(0).astype('int')
df['alignL_count'] = pd.to_numeric(df['alignL_count'], errors='coerce').fillna(0).astype('int')
df['fracMatchIndel_max'] = pd.to_numeric(df['fracMatchIndel_max'], errors='coerce').fillna(0.0).astype('float')

df = df.sample(frac=0.01)
train_cols = df.columns[0:7]
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

# INDELs


## INDELs with segdups and homopolymers

``` bash
bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels.bed -b GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed -b GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed

cat HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed | cut -f7,8,9,10,14,19 > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers.bed

cut -f1,2,3 HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_DP_VAF_indel_length.bed
cut -f4 HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_label.bed
cut -f5,6 HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_CG_AT_lens.bed
paste HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_DP_VAF_indel_length.bed HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_CG_AT_lens.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_DP_VAF_indel_length_CG_AT_lens_temp.bed
paste HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_DP_VAF_indel_length_CG_AT_lens_temp.bed HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_label.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_homopolymers_DP_VAF_indel_length_CG_AT_lens_label.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_indels_annotated_with_GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size_GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,10,14,19,25,26,33 > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers.bed

cut -f1,2,3 HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length.bed
cut -f4 HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_label.bed
cut -f5,6,7,8,9 HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_CG_AT_lens.bed
paste HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length.bed HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_CG_AT_lens.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length_CG_AT_lens_temp.bed
paste HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length_CG_AT_lens_temp.bed HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_label.bed > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length_CG_AT_lens_label.bed
```

### Pickling examples INDELs
``` python
import pickle 

from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

df = pd.read_csv("HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length_CG_AT_lens_label.bed", sep="\t")
df.columns = ["DP", "VAF", "indel_length", "CGhomopolgt3_len", "AThomopolgt3_len", "alignL_max", "alignL_count", "fracMatchIndel_max", "label"]

# adapted from #3 at https://towardsdatascience.com/converting-data-to-a-numeric-type-in-pandas-db9415caab0b
df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype('int')
df['VAF'] = pd.to_numeric(df['VAF'], errors='coerce').fillna(0.0).astype('float')
df['indel_length'] = pd.to_numeric(df['indel_length'], errors='coerce').fillna(0).astype('int')
df['CGhomopolgt3_len'] = pd.to_numeric(df['CGhomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['AThomopolgt3_len'] = pd.to_numeric(df['AThomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['alignL_max'] = pd.to_numeric(df['alignL_max'], errors='coerce').fillna(0).astype('int')
df['alignL_count'] = pd.to_numeric(df['alignL_count'], errors='coerce').fillna(0).astype('int')
df['fracMatchIndel_max'] = pd.to_numeric(df['fracMatchIndel_max'], errors='coerce').fillna(0.0).astype('float')

df = df.sample(frac=0.01)
train_cols = df.columns[0:8]
label = df.columns[-1]
X = df[train_cols]
y = df[label]

seed = 1

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=seed)

ebm = ExplainableBoostingClassifier(random_state=seed)
ebm.fit(X_train, y_train)

with open('one_percent_ebm_indels_data.pickle', 'wb') as f:
    pickle.dump(ebm, f, pickle.HIGHEST_PROTOCOL)

import pickle 
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show


with open('one_percent_ebm_indels_data.pickle', 'rb') as f:
     ebm = pickle.load(f)

ebm_global = ebm.explain_global()
show(ebm_global)
```


### Pickling examples SNPs
``` python
import pickle 
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

df = pd.read_csv("HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups_and_homopolymers_DP_VAF_CG_AT_lens_alignLMax_alignLCount_fracMatchIndelMax_label.bed", sep="\t")
df.columns = ["DP", "VAF", "CGhomopolgt3_len", "AThomopolgt3_len", "alignL_max", "alignL_count", "fracMatchIndel_max", "label"]

# adapted from #3 at https://towardsdatascience.com/converting-data-to-a-numeric-type-in-pandas-db9415caab0b
df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype('int')
df['VAF'] = pd.to_numeric(df['VAF'], errors='coerce').fillna(0.0).astype('float')
df['CGhomopolgt3_len'] = pd.to_numeric(df['CGhomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['AThomopolgt3_len'] = pd.to_numeric(df['AThomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['alignL_max'] = pd.to_numeric(df['alignL_max'], errors='coerce').fillna(0).astype('int')
df['alignL_count'] = pd.to_numeric(df['alignL_count'], errors='coerce').fillna(0).astype('int')
df['fracMatchIndel_max'] = pd.to_numeric(df['fracMatchIndel_max'], errors='coerce').fillna(0.0).astype('float')

df = df.sample(frac=0.01)
train_cols = df.columns[0:7]
label = df.columns[-1]
X = df[train_cols]
y = df[label]

seed = 1

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=seed)

ebm = ExplainableBoostingClassifier(random_state=seed)
ebm.fit(X_train, y_train)

with open('one_percent_ebm_snps_data.pickle', 'wb') as f:
    pickle.dump(ebm, f, pickle.HIGHEST_PROTOCOL)


import pickle 
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show


with open('one_percent_ebm_snps_data.pickle', 'rb') as f:
     ebm = pickle.load(f)

ebm_global = ebm.explain_global()
show(ebm_global)
```