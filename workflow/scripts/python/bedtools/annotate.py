from functools import reduce
from typing import List
import numpy as np
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import read_tsv, write_tsv
from common.cli import setup_logging
import common.config as cfg

logger = setup_logging(snakemake.log[0])  # type: ignore


def left_outer_intersect(left, path):
    logger.info("Adding annotations from %s", path)

    # Use bedtools to perform left-outer join of two bed/tsv files. Since
    # bedtools will join all columns from the two input files, keep track of the
    # width of the left input file so that the first three columns of the right
    # input (chr, chrStart, chrEnd, which are redundant) can be dropped.
    left_cols = left.columns.tolist()
    left_width = len(left_cols)
    right = read_tsv(path)
    # ASSUME the first three columns are the bed index columns
    right_cols = ["_" + c if i < 3 else c for i, c in enumerate(right.columns.tolist())]
    right_bed = bt.from_dataframe(right)
    # prevent weird type errors when converted back to dataframe from bed
    dtypes = {right_cols[0]: str}
    # convert "." to NaN since "." is a string/object which will make pandas run
    # slower than an actual panda
    na_vals = {c: "." for c in left_cols + right_cols[3:]}
    new_df = (
        bt.from_dataframe(left)
        .intersect(right_bed, loj=True)
        .to_dataframe(names=left_cols + right_cols, na_values=na_vals, dtype=dtypes)
    )
    # Bedtools intersect will use -1 for NULL in the case of numeric columns. I
    # suppose this makes sense since any "real" bed columns (according to the
    # "spec") will always be positive integers or strings. Since -1 might be a
    # real value and not a missing one in my case, use the chr field to figure
    # out if a row is "missing" and fill NaNs accordingly
    new_cols = new_df.columns[left_width:]
    new_pky = new_cols[:3]
    new_chr = new_pky[0]
    new_data_cols = new_cols[3:]
    new_df.loc[:, new_data_cols] = new_df[new_data_cols].where(
        new_df[new_chr] != ".", np.nan
    )

    logger.info("Annotations added: %s\n", ", ".join(new_data_cols))

    return new_df.drop(columns=new_pky)


def intersect_tsvs(config: cfg.StratoMod, ifile: str, ofile: str, tsv_paths: List[str]):
    target_df = read_tsv(ifile)
    new_df = reduce(left_outer_intersect, tsv_paths, target_df)
    new_df.insert(
        loc=0,
        column=config.feature_meta.raw_index,
        value=new_df.index,
    )
    write_tsv(ofile, new_df)


def main(smk, config: cfg.StratoMod) -> None:
    tsvs = smk.input.annotations
    vcf = smk.input.variants[0]
    logger.info("Adding annotations to %s\n", vcf)
    intersect_tsvs(config, vcf, smk.output[0], tsvs)


main(snakemake, snakemake.config)  # type: ignore
