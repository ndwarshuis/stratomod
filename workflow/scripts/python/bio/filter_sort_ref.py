import re
from typing import Any
import subprocess as sp
import common.config as cfg
from Bio import bgzf  # type: ignore
from common.io import setup_logging

logger = setup_logging(snakemake.log[0])  # type: ignore


def stream_fasta(ipath: str, chr_names: list[str]) -> sp.Popen[bytes]:
    return sp.Popen(
        ["samtools", "faidx", ipath, *chr_names],
        stdout=sp.PIPE,
        stderr=sp.PIPE,
    )


def stream_sdf(ipath: str, chr_names: list[str]) -> sp.Popen[bytes]:
    return sp.Popen(
        [
            *["rtg", "sdf2fasta", "--no-gzip", "--line-length=70"],
            *["--input", ipath],
            *["--output", "-"],
            *["--names", *chr_names],
        ],
        stdout=sp.PIPE,
        stderr=sp.PIPE,
    )


def main(smk: Any, sconf: cfg.StratoMod) -> None:
    rsk = cfg.RefsetKey(smk.wildcards["refset_key"])
    cs = sconf.refsetkey_to_chr_indices(rsk)
    prefix = sconf.refsetkey_to_ref(rsk).sdf.chr_prefix

    chr_mapper = {c.chr_name_full(prefix): c.value for c in cs}
    chr_names = [*chr_mapper]

    # Read from a fasta or sdf depending on what we were given; in either
    # case, read only the chromosomes we want in sorted order and return a
    # fasta text stream
    def choose_input(i: Any) -> sp.Popen[bytes]:
        try:
            return stream_fasta(i.fasta[0], chr_names)
        except AttributeError:
            try:
                return stream_sdf(i.sdf[0], chr_names)
            except AttributeError:
                assert False, "unknown input key, this should not happen"

    p = choose_input(smk.input)

    if p.stdout is not None:
        # Stream the fasta and replace the chromosome names in the header with
        # its integer index
        with bgzf.open(smk.output[0], "w") as f:
            for i in p.stdout:
                if i.startswith(b">"):
                    m = re.match(">([^ \n]+)", i.decode())
                    if m is None:
                        logger.error("could get chrom name from FASTA header")
                        exit(1)
                    try:
                        f.write(f">{chr_mapper[m[1]]}\n")
                    except KeyError:
                        assert False, (
                            "could not convert '%s' to index, this should not happen"
                            % m[1]
                        )
                else:
                    f.write(i)
    else:
        assert False, "stdout not a pipe, this should not happen"

    p.wait()

    if p.returncode != 0:
        logger.error(p.stderr)
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
