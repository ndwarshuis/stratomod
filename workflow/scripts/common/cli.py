import argparse


def add_input_arg(desc, parser):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help=(
            "the path from which to read {}; " "if not given, read from stdin"
        ).format(desc),
    )


def add_output_arg(desc, parser):
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help=("the path to which to write {};" "if not given, write to stdout").format(
            desc
        ),
    )


def add_config_arg(desc, parser):
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="{} as a JSON string".format(desc),
    )


def make_io_parser(desc, idesc, odesc):
    parser = argparse.ArgumentParser(description=desc)
    add_output_arg(odesc, parser)
    add_input_arg(idesc, parser)
    return parser
