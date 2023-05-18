from typing import Any, TextIO
from common.config import StratoMod, RefsetKey
from common.io import with_gzip_maybe


def filter_file(smk: Any, config: StratoMod, fi: TextIO, fo: TextIO) -> None:
    rsk = RefsetKey(smk.wildcards["refset_key"])
    chr_prefix = smk.params.chr_prefix
    cs = config.refsetkey_to_chr_indices(rsk)

    chr_mapper = {c.chr_name_full(chr_prefix): c.value for c in cs}

    for ln in fi:
        if ln.startswith("#"):
            fo.write(ln)
        else:
            ls = ln.rstrip().split("\t")
            try:
                ls[0] = str(chr_mapper[ls[0]])
                fo.write("\t".join(ls) + "\n")
            except KeyError:
                pass


def main(smk: Any, config: StratoMod) -> None:
    with_gzip_maybe(
        lambda i, o: filter_file(smk, config, i, o),
        str(smk.input[0]),
        str(smk.output[0]),
    )


main(snakemake, snakemake.config)  # type: ignore
