#! /usr/bin/env python3

import sys
import yaml
from workflow.scripts.python.common.config import StratoMod

path = sys.argv[1]
modkey = sys.argv[2]

with open(path, "r") as f:
    config = StratoMod.parse_obj(yaml.safe_load(f))
    fs = sorted(config.modelkey_to_features(modkey).items(), key=lambda x: x[0])
    for n, d in fs:
        print(f"{n}: {d}\n\n")
