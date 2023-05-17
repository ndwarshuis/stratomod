# Features Overview

The following is a list and short description of each feature in the annotated
dateframe (eg that which should be output to
`results/annotated_input/*/{SNP,INDEL}.tsv`).

## Primary Key Fields

Each row in the dataframe is one variant. Each variant is uniquely identified by
the chromosome on which it is located (`CHROM`), and the start/end positions on
this chromosome (`POS` and `POS+length(REF)`).

`CHROM` is a string like `chrN` where `N` is a number 1-22 or `X` or `Y`. The
two position fields are zero-indexed integers denoting the starting offset in
the chromosome.

Each feature described is also denoted by the same three fields. In order for a
feature to be 'part of a variant' the region of the feature must overlap with
that of the variant (like a SQL join but where the join condition is
'overlapping region' and not 'equality.'

## Label

The label field is simply called `label`. This has one of 4 values:

* tp: variant is in both truth and query
* fp: variant is in only the query
* fn: variant is in only the truth
* tn: variant is in neither the truth or query (which will only happen when we
  count filtered variants as not part of the query, in which case the variant
  will still get a row if it is not also in the truth)

Since the EBMs are used a binary classifer, these labels need to be
filtered/mapped to 0 or 1 before being fed to the EBM.

## Genomic context features

These vary depending on the model being run. To see a list of valid features
for a given model key, run the included script (with the stratomod conda env
activated):

```
dump_features path/to/config.yml <model_key>
```

In addition to being present in the annotated data frame, any one of these
features may be specified in the `features` block in a given model.
