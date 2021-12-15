# simple repeat regions

Generate a bed file with non-overlapping regions with summary statistics for
simple repeats. See
[here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema).

## dependencies

- python == 3.9
- pandas
- pybedtools

## obtaining files

### simple repeat table

``` bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
      -A -P 3306 -D hg38 -B \
      -e 'select * from simpleRepeat;' \
      > simple_repeat.tsv
```

### genome file

This is necessary to add slop

``` bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
      -A -P 3306 -D hg38 -N -B \
      -e 'select chrom,size from chromInfo;' \
      > genome.txt
```

## generating merged bed file

See/run this:

``` bash
python3 merge_repeats.py
```
