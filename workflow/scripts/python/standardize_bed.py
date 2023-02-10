import gzip
import io
import common.config as cfg


def filter_file(smk, fi: io.TextIOWrapper) -> None:
    input_key = smk.wildcards["input_key"]
    chr_prefix = cfg.inputkey_to_chr_prefix(smk.config, input_key)
    chr_filter = cfg.inputkey_to_chr_filter(smk.config, input_key)
    fs = tuple(["#", *[f"{c}\t" for c in chr_filter]])

    def tolines(f: io.TextIOWrapper):
        f.writelines(
            (y for x in fi if (y := x.removeprefix(chr_prefix)).startswith(fs)),
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
