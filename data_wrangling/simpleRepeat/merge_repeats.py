import pandas as pd
from itertools import repeat
from pybedtools import BedTool

srs_df = pd.read_csv("simple_repeat.tsv", sep="\t").drop(columns=["bin"])
srs_columns = srs_df.columns.tolist()

srs_bed = BedTool.from_dataframe(srs_df)

col_indices = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]


def repeat_opt(o):
    return list(repeat(o, len(col_indices)))


# just use one column for count since all columns will produce the same number
opts = ["count"] + repeat_opt("max") + repeat_opt("min") + repeat_opt("median")
cols = [5] + col_indices + col_indices + col_indices

headers = ["chr", "chromStart", "chromEnd", "count"] + [
    "{}_{}".format(srs_columns[c - 1], o) for c, o in list(zip(cols, opts))[1:]
]

merged = srs_bed.merge(c=cols, o=opts).slop(b=5, g="genome.txt")

merged_df = merged.to_dataframe(names=headers)
# use the original dataframe to get the region length since we added slop to the
# merged version
merged_df["region_length"] = srs_df["chromEnd"] - srs_df["chromStart"]

merged_df.to_csv("simple_repeat_merged.bed", sep="\t", index=False)
