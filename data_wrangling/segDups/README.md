# segmental duplications

Annotate with genomicSuperDups.txt downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

## dependencies

- bedtools => v2.27.1

``` bash
cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 4 -o min,max,count,mean > genomicSuperDups_alignL_stats.bed

cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 5 -o min,max,count,mean > genomicSuperDups_fracMatchIndel_stats.bed
```

## generating bed file with annotations for alignL and fracMatchIndel

```
# SNPs
bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_snps_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,13,14,15,16,21,22,23,24 > HG002_DV_Illumina_df_ebm_snps_annotated_with_segdups.bed

# INDELs
bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels.bed -b genomicSuperDups_alignL_stats.bed > HG002_DV_Illumina_df_ebm_indels_genomicSuperDups_alignL_stats.bed

bedtools intersect -wao -a HG002_DV_Illumina_df_ebm_indels_genomicSuperDups_alignL_stats.bed -b genomicSuperDups_fracMatchIndel_stats.bed > HG002_DV_Illumina_df_ebm_indels_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed

cat HG002_DV_Illumina_df_ebm_indels_genomicSuperDups_alignL_stats_genomicSuperDups_fracMatchIndel_stats.bed | cut -f7,8,9,13,14,15,16,21,22,23,24 > HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups.bed
```