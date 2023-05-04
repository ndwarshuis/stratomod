import gzip
import hashlib
from pathlib import Path
from logging import Logger


def get_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with gzip.open(p, "r") as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False


# def is_bgzip(p: Path) -> bool:
#     # since bgzip is in blocks (vs gzip), determine if in bgzip by
#     # attempting to seek first block
#     with open(p, "rb") as f:
#         try:
#             next(bgzf.BgzfBlocks(f), None)
#             return True
#         except ValueError:
#             return False


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
