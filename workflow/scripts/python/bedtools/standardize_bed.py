import gzip
import io
import re
from typing import Callable
from common.config import StratoMod, ChrIndex

# I could use pandas for all this, but vcfeval will complain if I strip out the
# headers (which pandas will do). Imperative loop it is...


def filter_file(smk, config: StratoMod, fi: io.TextIOWrapper) -> None:
    chr_prefix = smk.params.chr_prefix
    chr_indices = config.refsetkey_to_chr_indices(smk.wildcards["refset_key"])
    fs = tuple(["#", *[f"{i.chr_name}\t" for i in chr_indices]])

    def make_sub(i: ChrIndex) -> Callable[[str], str]:
        pat = re.compile(f"^{i.chr_name}")
        return lambda s: pat.sub(str(i.value), s)

    subX = make_sub(ChrIndex.CHRX)
    subY = make_sub(ChrIndex.CHRY)

    # pre 3.9 :(
    def remove_prefix(s: str) -> str:
        if s.startswith(chr_prefix):
            return s[len(chr_prefix) :]  # noqa: ignore=E203
        return s

    def tolines(f: io.TextIOWrapper):
        f.writelines(
            (subY(subX(y)) for x in fi if (y := remove_prefix(x)).startswith(fs)),
        )

    if str(smk.output[0]).endswith(".gz"):
        with gzip.open(smk.output[0], "wt") as fo:
            tolines(fo)
    else:
        with open(smk.output[0], "wt") as fo:
            tolines(fo)


def main(smk, config: StratoMod) -> None:
    if str(smk.input[0]).endswith(".gz"):
        with gzip.open(smk.input[0], "rt") as f:
            filter_file(smk, config, f)
    else:
        with open(smk.input[0], "rt") as f:
            filter_file(smk, config, f)


main(snakemake, snakemake.config)  # type: ignore
