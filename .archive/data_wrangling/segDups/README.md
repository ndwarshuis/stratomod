# segmental duplications

Annotate with genomicSuperDups.txt downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

## dependencies

- bedtools => v2.27.1

``` bash
cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 4 -o min,max,count,mean > genomicSuperDups_alignL_stats.bed

cat genomicSuperDups.txt | cut -f2-4,19,28 | mergeBed -i stdin -c 5 -o min,max,count,mean > genomicSuperDups_fracMatchIndel_stats.bed
```