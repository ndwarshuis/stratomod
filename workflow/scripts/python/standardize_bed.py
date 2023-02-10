import gzip
import io
import common.config as cfg


def filter_file(fi: io.TextIOWrapper) -> None:
    chr_prefix = cfg.inputkey_to_chr_prefix(
        snakemake.config,
        snakemake.wildcards["input_key"],
    )
    chr_filter = cfg.inputkey_to_chr_filter(
        snakemake.config,
        snakemake.wildcards["input_key"],
    )
    contig = "##contig=<ID="
    contig_prefix = contig + chr_prefix
    fs = tuple(["#", *[f"{c}\t" for c in chr_filter]])
    with open(snakemake.output[0], "w") as fo:
        fo.writelines(
            x.replace(contig_prefix, contig) if x.startswith("#") else x
            for x in filter(
                lambda x: x.startswith(fs),
                (x.removeprefix(chr_prefix) for x in fi.readlines()),
            )
        )


def main() -> None:
    if snakemake.params.unzip is True:
        with gzip.open(snakemake.input[0], "rt") as f:
            filter_file(f)
    else:
        with open(snakemake.input[0], "rt") as f:
            filter_file(f)


main()
