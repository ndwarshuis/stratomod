import gzip
import io
import common.config as cfg


def filter_file(smk, fi: io.TextIOWrapper) -> None:
    chr_prefix = cfg.inputkey_to_chr_prefix(
        smk.config,
        smk.wildcards["input_key"],
    )
    chr_filter = cfg.inputkey_to_chr_filter(
        smk.config,
        smk.wildcards["input_key"],
    )
    # contig = "##contig=<ID="
    # contig_prefix = contig + chr_prefix
    fs = tuple(["#", *[f"{c}\t" for c in chr_filter]])

    def tolines(f: io.TextIOWrapper):
        f.writelines(
            # x.replace(contig_prefix, contig) if x.startswith("#") else x
            x
            for x in filter(
                lambda x: x.startswith(fs),
                (x.removeprefix(chr_prefix) for x in fi.readlines()),
            )
        )

    if smk.params.gzip_out is True:
        with gzip.open(smk.output[0], "wt") as fo:
            tolines(fo)
    else:
        with open(smk.output[0], "wt") as fo:
            tolines(fo)


def main(smk) -> None:
    if smk.params.gzip_in is True:
        with gzip.open(smk.input[0], "rt") as f:
            filter_file(smk, f)
    else:
        with open(smk.input[0], "rt") as f:
            filter_file(smk, f)


main(snakemake)  # type: ignore
