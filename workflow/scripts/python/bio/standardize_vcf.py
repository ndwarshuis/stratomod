import gzip
from typing import Any, cast, TextIO
from common.config import StratoMod, RefsetKey, VCFFile
from Bio import bgzf  # type: ignore


def fix_DV_refcall(filter_col: str, sample_col: str) -> str:
    return (
        sample_col.replace("./.", "0/1").replace("0/0", "0/1")
        if filter_col == "RefCall"
        else sample_col
    )


def strip_format_fields(
    fields: set[str],
    format_col: str,
    sample_col: str,
) -> tuple[str, str]:
    f, s = zip(
        *[
            (f, s)
            for f, s in zip(format_col.split(":"), sample_col.split(":"))
            if f not in fields
        ]
    )
    return (":".join(f), ":".join(s))


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
            ls = ln.rstrip().split("\t")[:10]
            # CHROM = 0
            # POS = 1
            # ID = 2
            # REF = 3
            # ALT = 4
            # QUAL = 5
            # FILTER = 6
            # INFO = 7
            # FORMAT = 8
            # SAMPLE = 9
            try:
                ls[0] = str(chr_mapper[ls[0]])
                if vcf.corrections.fix_refcall_gt:
                    ls[9] = fix_DV_refcall(ls[6], ls[9])
                if len(vcf.corrections.strip_format_fields) > 0:
                    ls[8], ls[9] = strip_format_fields(
                        vcf.corrections.strip_format_fields,
                        ls[8],
                        ls[9],
                    )
                fo.write("\t".join(ls) + "\n")
            except KeyError:
                pass


def main(smk: Any, config: StratoMod) -> None:
    i = str(smk.input[0])
    o = str(smk.output[0])
    # bgzf only understands latin1, so read everything as such
    with gzip.open(i, "rt", encoding="latin1") if i.endswith(".gz") else open(
        i, "rt", encoding="latin1"
    ) as fi:
        with bgzf.open(o, "wt") if o.endswith(".gz") else open(o, "wt") as fo:
            filter_file(smk, config, fi, fo)


main(snakemake, snakemake.config)  # type: ignore
