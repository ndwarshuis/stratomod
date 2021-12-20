import pandas as pd
from pybedtools import BedTool

rmsk_df = pd.read_csv("rmsk.tsv", sep="\t")

classcol = "repClass"


def merge_class(df, classname):
    dropped = df[df[classcol] == classname].drop(columns=[classcol])
    merged = BedTool.from_dataframe(dropped).merge().to_dataframe()
    merged["length"] = merged["end"] - merged["start"]
    merged[classcol] = classname
    return merged


rmsk_merged = pd.concat(
    [
        merge_class(rmsk_df, "SINE"),
        merge_class(rmsk_df, "LINE"),
        merge_class(rmsk_df, "LTR"),
        merge_class(rmsk_df, "Satellite"),
    ],
    axis=0,
)

rmsk_merged.to_csv("rmsk_merged.bed", sep="\t", index=False)
