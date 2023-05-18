import gzip
import hashlib
from typing import Callable, TextIO, TypeVar
from pathlib import Path
from logging import Logger

X = TypeVar("X")


def with_gzip_maybe(f: Callable[[TextIO, TextIO], X], i: str, o: str) -> X:
    hi = (
        gzip.open(i, "rt", encoding="latin1")
        if i.endswith(".gz")
        else open(i, "rt", encoding="latin1")
    )
    ho = gzip.open(o, "wt") if o.endswith(".gz") else open(o, "wt")
    with hi as fi, ho as fo:
        return f(fi, fo)


def get_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def get_md5_dir(path: Path) -> str:
    """Get the "MD5" of a directory.

    This is equivalent to the following bash command:

    find <dir> -type f -exec md5sum {} \\; | \\
        grep -v snakemake_timestamp | \\
        LC_COLLATE=C sort -k 2 | \\
        cut -d' ' -f1 | \\
        tr -d '\n' | \\
        md5sum

    Will only get the md5 of the first layer of files (eg not recursive). This
    is appropriate for SDF files which have this property.
    """
    h = hashlib.md5()
    ps = sorted(
        [
            p
            for p in path.iterdir()
            if p.is_file() and "snakemake_timestamp" not in p.name
        ]
    )
    for p in ps:
        h.update(get_md5(p).encode())
    return h.hexdigest()


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with gzip.open(p, "r") as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False


# set up basic logger that prints to both console and a file (the log directive
# from snakemake) and captures warnings so those don't go unnoticed
def setup_logging(path: str, console: bool = False) -> Logger:
    import logging

    logging.basicConfig(filename=path, level=logging.INFO)
    logging.captureWarnings(True)
    logger = logging.getLogger()
    if console:
        logger.addHandler(logging.StreamHandler())
    return logger
