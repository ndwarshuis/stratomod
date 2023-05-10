from pathlib import Path
import subprocess as sp
from typing import Callable, Any, cast
from typing_extensions import assert_never
from tempfile import NamedTemporaryFile as Tmp
from common.io import is_gzip, setup_logging, get_md5, get_md5_dir
from common.bed import is_bgzip
import common.config as cfg


GUNZIP = ["gunzip", "-c"]
BGZIP = ["bgzip", "-c"]
CURL = ["curl", "-Ss", "-L", "-q"]

log = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, params: Any) -> None:
    src = cast(params.src, cfg.FileSrc)
    opath = Path(smk.output[0])
    is_fasta = smk.params.is_fasta

    if isinstance(src, cfg.LocalSrc):
        # ASSUME this is in the format we indicate (TODO be more paranoid)
        opath.symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HTTPSrc):
        curlcmd = [*CURL, src.url]

        if is_fasta:
            # to test the format of downloaded files, sample the first 65000 bytes
            # (which should be enough to get one block of a bgzip file, which will
            # allow us to test for it)
            curltestcmd = [*CURL, "-r", "0-65000", src.url]

            with open(opath, "wb") as f, Tmp() as tf:

                def curl() -> None:
                    sp.Popen(curlcmd, stdout=f).wait()

                def curl_test(testfun: Callable[[Path], bool]) -> bool:
                    sp.Popen(curltestcmd, stdout=tf).wait()
                    return testfun(Path(tf.name))

                def curl_gzip(cmd: list[str]) -> None:
                    p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                    p2 = sp.Popen(cmd, stdin=p1.stdout, stdout=f)
                    p2.wait()

                if curl_test(is_bgzip):
                    curl()
                elif curl_test(is_gzip):
                    p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                    p2 = sp.Popen(GUNZIP, stdin=p1.stdout, stdout=sp.PIPE)
                    p3 = sp.Popen(BGZIP, stdin=p2.stdout, stdout=f)
                    p3.wait()
                else:
                    curl_gzip(BGZIP)

        else:
            tarcmd = [
                *["bsdtar", "-xf", "-"],
                *["--directory", str(opath)],
                "--strip-component=1",
            ]

            opath.mkdir(parents=True)

            p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
            p2 = sp.Popen(tarcmd, stdin=p1.stdout)
            p2.wait()

    else:
        assert_never(src)

    if src.md5 is not None:
        actual = get_md5(opath) if is_fasta else get_md5_dir(opath)
        if actual != src.md5:
            log.error("md5s don't match; wanted %s, actual %s", src.md5, actual)
            exit(1)


main(snakemake.output, snakemake.params)  # type: ignore
