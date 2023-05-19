#! /usr/bin/env python3

import sys
import yaml
from textwrap import fill
from workflow.scripts.python.common.config import StratoMod

path = sys.argv[1]
modkey = sys.argv[2]

with open(path, "r") as f:
    config = StratoMod.parse_obj(yaml.safe_load(f))
    fs = sorted(config.modelkey_to_features(modkey).items(), key=lambda x: x[0])
    for prefix, (desc0, ss) in fs:
        print(f"Feature Prefix: {prefix}")
        print(
            fill(
                f"Description: {desc0}",
                width=70,
                initial_indent="  ",
                subsequent_indent="  ",
                break_long_words=False,
            )
        )
        print()
        print(f"  Features:")
        for feat, desc1 in ss.items():
            print(
                fill(
                    f"{feat}: {desc1}",
                    width=70,
                    initial_indent="  - ",
                    subsequent_indent="    ",
                )
            )
            print()
