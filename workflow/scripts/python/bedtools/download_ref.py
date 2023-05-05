from pathlib import Path
import subprocess as sp
from typing_extensions import assert_never
import common.config as cfg

# from typing import Callable
# from tempfile import NamedTemporaryFile as Tmp
from common.io import get_md5  # , is_gzip, is_bgzip


# log = setup_logging(snakemake.log[0])  # type: ignore

GUNZIP = ["gunzip", "-c"]
BGZIP = ["bgzip", "-c"]
CURL = ["curl", "-Ss", "-L", "-q"]


def main(opath: Path, src: cfg.FileSrc) -> None:
    if isinstance(src, cfg.LocalSrc):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        opath.symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HTTPSrc):
        curlcmd = [*CURL, src.url]
        tarcmd = [
            "bsdtar",
            "-xf",
            "-",
            "--directory",
            str(opath),
            "--strip-component=1",
        ]

        opath.mkdir(parents=True)

        p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
        p2 = sp.Popen(tarcmd, stdin=p1.stdout)
        p2.wait()

        # to test the format of downloaded files, sample the first 65000 bytes
        # (which should be enough to get one block of a bgzip file, which will
        # allow us to test for it)
        # curltestcmd = [*CURL, "-r", "0-65000", src.url]

        # with Tmp() as tf:

        #     def curl_test(testfun: Callable[[Path], bool]) -> bool:
        #         sp.Popen(curltestcmd, stdout=tf).wait()
        #         return testfun(Path(tf.name))

        #     if curl_test(is_bgzip):
        #         untar(p1)
        #     elif curl_test(is_gzip):
        #         p2 = sp.Popen(GUNZIP, stdin=p1.stdout, stdout=sp.PIPE)
        #         untar(sp.Popen(BGZIP, stdin=p2.stdout, stdout=sp.PIPE))
        #     else:
        #         untar(sp.Popen(BGZIP, stdin=p1.stdout, stdout=sp.PIPE))

    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := get_md5(opath)):
        # log.error("md5s don't match; wanted %s, actual %s", src.md5, actual)
        # TODO actually log this
        print("md5s don't match; wanted %s, actual %s" % (src.md5, actual))
        exit(1)


main(Path(snakemake.output[0]), snakemake.params.src)  # type: ignore
