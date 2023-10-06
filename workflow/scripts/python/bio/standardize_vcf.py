import re
from pathlib import Path
import subprocess as sp
from typing import Any, cast, IO, NamedTuple
from common.config import StratoMod, RefsetKey, VCFFile
from common.bed import with_bgzip_maybe
from common.io import setup_logging

log = setup_logging(snakemake.log["filtered"])  # type: ignore


class Env(NamedTuple):
    refsetkey: RefsetKey
    vcf: VCFFile
    config: StratoMod
    split_log: Path


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


def filter_file(env: Env, fi: IO[str], fo: IO[str]) -> None:
    chr_prefix = env.vcf.chr_prefix
    cs = env.config.refsetkey_to_chr_indices(env.refsetkey)

    chr_mapper = {c.chr_name_full(chr_prefix): c.value for c in cs}

    for ln in fi:
        if ln.startswith("#"):
            if cmatch := re.match("##contig=<ID=([^,]+),length=(\\d+)>", ln):
                try:
                    chrom_id = chr_mapper[cmatch[1]]
                    fo.write(f"##contig=<ID={chrom_id},length={cmatch[2]}>\n")
                except KeyError:
                    pass
            else:
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
                if env.vcf.corrections.fix_refcall_gt:
                    ls[9] = fix_DV_refcall(ls[6], ls[9])
                if len(env.vcf.corrections.strip_format_fields) > 0:
                    ls[8], ls[9] = strip_format_fields(
                        env.vcf.corrections.strip_format_fields,
                        ls[8],
                        ls[9],
                    )
                fo.write("\t".join(ls) + "\n")
            except KeyError:
                pass


def split_biallelics_and_filter_file(env: Env, fi: IO[str], fo: IO[str]) -> None:
    with open(env.split_log, "wt") as lo:
        cmd = ["bcftools", "norm", "--multiallelics", "-"]
        p = sp.Popen(cmd, stdin=sp.PIPE, stdout=fo, stderr=lo, text=True)
        if p.stdin is not None:
            filter_file(env, fi, p.stdin)
            # manually close stdin to let the process exit; I could do this
            # also probably by letting the context managers scope out but then
            # don't get a nice containerized function
            p.stdin.close()
            p.wait()
            if p.returncode != 0:
                log.error("error when splitting bialellics")
        else:
            assert False, "this should not happen"


def main(smk: Any, config: StratoMod) -> None:
    env = Env(
        refsetkey=RefsetKey(smk.wildcards["refset_key"]),
        vcf=cast(VCFFile, smk.params.vcf),
        config=config,
        split_log=smk.log["split"],
    )

    if env.vcf.split_biallelics:
        action = split_biallelics_and_filter_file
    else:
        action = filter_file

    with_bgzip_maybe(
        lambda i, o: action(env, i, o),
        str(smk.input[0]),
        str(smk.output[0]),
    )


main(snakemake, snakemake.config)  # type: ignore
