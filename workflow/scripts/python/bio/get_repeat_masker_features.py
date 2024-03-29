import pandas as pd
from pathlib import Path
from typing import Optional, Any
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.tsv import write_tsv
from common.io import setup_logging
from common.bed import read_bed

# The repeat masker database is documented here:
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

logger = setup_logging(snakemake.log[0])  # type: ignore

# both of these columns are temporary and used to make processing easier
CLASSCOL = "_repClass"
FAMCOL = "_repFamily"


def main(smk: Any, config: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    rk = config.refsetkey_to_refkey(rsk)
    src = config.references[rk].feature_data.repeat_masker
    cs = config.refsetkey_to_chr_indices(rsk)

    def read_rmsk_df(path: Path) -> pd.DataFrame:
        cols = {11: CLASSCOL, 12: FAMCOL}
        return read_bed(path, src.params, cols, cs)

    def merge_and_write_group(
        df: pd.DataFrame,
        groupcol: str,
        clsname: str,
        famname: Optional[str] = None,
    ) -> None:
        groupname = clsname if famname is None else famname
        dropped = df[df[groupcol] == groupname].drop(columns=[groupcol])
        merged = bt.from_dataframe(dropped).merge().to_dataframe(names=cfg.BED_COLS)
        col = config.feature_definitions.repeat_masker.fmt_name(src, clsname, famname)
        merged[col] = merged[cfg.BED_END] - merged[cfg.BED_START]
        write_tsv(smk.output[0], merged, header=True)

    cls = smk.wildcards.rmsk_class
    df = read_rmsk_df(smk.input[0])

    try:
        fam = smk.wildcards.rmsk_family
        logger.info("Filtering/merging rmsk family %s/class %s", fam, cls)
        merge_and_write_group(df, FAMCOL, cls, fam)
    except AttributeError:
        logger.info("Filtering/merging rmsk class %s", cls)
        merge_and_write_group(df, CLASSCOL, cls)


main(snakemake, snakemake.config)  # type: ignore
