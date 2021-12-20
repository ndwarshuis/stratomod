# repeat masker regions

Generate a bed file with non-overlapping regions for SINE, LINE, Satellite, and
LTR repeats in the repeat masker database.
See
[here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema).

## dependencies

- python == 3.9
- pandas
- pybedtools

## obtaining files

### repeat masker table

``` bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
      -A -P 3306 -D hg38 -B \
      -e 'select genoName,genoStart,genoEnd,repClass from rmsk;' \
      > rmsk.tsv
```

## generating merged bed file

See/run this:

``` bash
python3 merge_rmsk.py
```
