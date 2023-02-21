import pandas as pd
from common.tsv import read_tsv, write_tsv
from common.bed import sort_bed_numerically

# use pandas here since it will more reliably account for headers
df = pd.concat([read_tsv(i, header=0) for i in snakemake.input])  # type: ignore
write_tsv(snakemake.output[0], sort_bed_numerically(df))  # type: ignore
