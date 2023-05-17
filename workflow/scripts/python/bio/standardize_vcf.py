import gzip
from typing import Any, TextIO, cast
from common.config import StratoMod, RefsetKey, VCFFile


def fix_DV_refcall(filter_col: str, sample_col: str) -> str:
    return (
        sample_col.replace("./.", "0/1").replace("0/0", "0/1")
        if filter_col == "RefCall"
        else sample_col
    )


def strip_IPS(format_col: str, sample_col: str) -> tuple[str, str]:
    if "IPS" in format_col:
        f, s = zip(
            *[
                (f, s)
                for f, s in zip(":".split(format_col), ":".split(sample_col))
                if f != "IPS"
            ]
        )
        return (":".join(f), ":".join(s))
    else:
        return (format_col, sample_col)


def filter_file(smk: Any, config: StratoMod, fi: TextIO, fo: TextIO) -> None:
    rsk = RefsetKey(smk.wildcards["refset_key"])
    vcf = cast(VCFFile, smk.params.vcf)
    chr_prefix = vcf.chr_prefix
    cs = config.refsetkey_to_chr_indices(rsk)

    chr_mapper = {c.chr_name_full(chr_prefix): c.value for c in cs}

    for ln in fi:
        if ln.startswith("#"):
            fo.write(ln)
        else:
            ls = "\t".split(ln)
            try:
                ls[0] = chr_mapper[ls[0]]
                if vcf.fix_refcall:
                    ls[9] = fix_DV_refcall(ls[8], ls[9])
                if vcf.strip_IPS:
                    ls[8], ls[9] = strip_IPS(ls[8], ls[9])
                fo.write("\t".join(ls) + "\n")
            except KeyError:
                pass


def main(smk: Any, config: StratoMod) -> None:
    i = str(smk.input[0])
    o = str(smk.output[0])
    with gzip.open(i, "rt") if i.endswith(".gz") else open(i, "rt") as fi:
        with gzip.open(o, "wt") if o.endswith(".gz") else open(o, "wt") as fo:
            filter_file(smk, config, fi, fo)


main(snakemake, snakemake.config)  # type: ignore
