import re
from typing import Any
import subprocess as sp
import common.config as cfg
from Bio import bgzf  # type: ignore


def stream_fasta(ipath: str, chr_names: list[str]) -> sp.Popen[bytes]:
    return sp.Popen(["samtools", "faidx", ipath, *chr_names], stdout=sp.PIPE)


def stream_sdf(ipath: str, chr_names: list[str]) -> sp.Popen[bytes]:
    return sp.Popen(
        [
            *["rtg", "sdf2fasta", "-Z", "--line_length=70"],
            *["-i", ipath],
            *["-o", "-"],
            *["-n", *chr_names],
        ],
        stdout=sp.PIPE,
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
            return stream_fasta(i.fasta, chr_names)
        except AttributeError:
            try:
                return stream_sdf(i.sdf, chr_names)
            except AttributeError:
                assert False, "unknown input key, this should not happen"

    p = choose_input(smk.input)

    if p.stdout is not None and p.returncode == 0:
        # Stream the fasta and replace the chromosome names in the header with
        # its integer index
        with bgzf.open(smk.output[0], "w") as f:
            for i in p.stdout:
                if i.startswith(b">"):
                    m = re.match(">[^ ]+", i.decode())
                    assert m is not None, "could get chrom name from FASTA header"
                    try:
                        f.write(f">{chr_mapper[m[1]]}")
                    except KeyError:
                        assert False, "could not convert '%s' to index" % m[1]
                else:
                    f.write(i)
    else:
        assert False  # TODO make useful


main(snakemake, snakemake.config)  # type: ignore
