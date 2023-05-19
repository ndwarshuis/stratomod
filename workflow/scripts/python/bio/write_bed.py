import common.config as cfg
from Bio import bgzf  # type: ignore


def main(opath: str, regions: list[cfg.BedRegion]) -> None:
    with bgzf.open(opath, "w") as f:
        for r in sorted(regions):
            f.write(r.fmt() + "\n")


main(snakemake.output[0], snakemake.params.regions)  # type: ignore
