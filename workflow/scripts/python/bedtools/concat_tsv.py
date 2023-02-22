import pandas as pd
from common.tsv import write_tsv
from common.bed import sort_bed_numerically

# use pandas here since it will more reliably account for headers
df = pd.concat([pd.read_table(i, header=0) for i in snakemake.input])  # type: ignore
write_tsv(snakemake.output[0], sort_bed_numerically(df))  # type: ignore
