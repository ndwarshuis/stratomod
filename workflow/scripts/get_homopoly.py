import pandas as pd
from pybedtools import BedTool as bt
from pybedtools import cleanup
from common.tsv import write_tsv, read_tsv
from common.bed import sort_bed_numerically

START_COL = "start"
END_COL = "end"
BED_COLS = ["chr", START_COL, END_COL]

BASE_COL = "base"
SIMPLE_REPEAT_BED_COLS = [*BED_COLS, BASE_COL]

IMP_LEN_COL = "imperfect_length"
IMP_BED_COLS = [*BED_COLS, IMP_LEN_COL]

SLOP = 5

# TODO this works but means we will sort a massive dataframe each time the
# script is called
def read_input(path):
    print(f"Reading dataframe from {path}")
    df = read_tsv(path, header=None, comment="#", names=SIMPLE_REPEAT_BED_COLS)
    return sort_bed_numerically(df)


def filter_base(df, base, genome):
    print(f"Filtering bed file for {base}s")
    _df = df[df[BASE_COL] == f"unit={base}"].drop(columns=[BASE_COL])
    ldf = len(_df)
    assert ldf > 0, f"Filtered bed file for {base} has no rows"
    print(f"Merging {ldf} rows for {base}s and adding {SLOP}bp slop")
    _df[IMP_LEN_COL] = _df[END_COL] - _df[START_COL]
    merged = (
        bt.from_dataframe(_df)
        .merge(d=1, c=[4], o=["sum"])
        .slop(b=SLOP, g=genome)
        .to_dataframe(names=IMP_BED_COLS)
    )
    # these files are huge; now that we have a dataframe, remove all the bed
    # files from tmpfs to prevent a run on downloadmoreram.com
    cleanup()
    return merged


def filter_bases(df, bases, genome):
    return [filter_base(df, b, genome) for b in bases]


def intersect_bases(dfs, bases):
    print(f"Concatenating and sorting merged beds for {bases}")
    sdf = sort_bed_numerically(pd.concat(dfs))

    print(f"Merging {bases} homopolymers")
    merged = (
        bt.from_dataframe(sdf)
        .merge(
            c=[4],
            o=["sum"],
        )
        .to_dataframe(names=IMP_BED_COLS)
    )

    print("Adding homopolymer length and perfect fraction")
    length_col = f"{bases}_homopolymer_length"
    frac_col = f"{bases}_homopolymer_perfect_frac"
    merged[length_col] = merged[END_COL] - merged[START_COL] - SLOP * 2
    merged[frac_col] = merged[IMP_LEN_COL] / merged[length_col]
    return merged.drop(columns=[IMP_LEN_COL])


def main():
    # ASSUME this file is already sorted
    simreps = read_input(snakemake.input["bed"][0])
    bases = snakemake.wildcards["bases"]
    base_dfs = filter_bases(simreps, bases, snakemake.input["genome"][0])
    merged_bases = intersect_bases(base_dfs, bases)
    write_tsv(snakemake.output[0], merged_bases, header=True)


main()
