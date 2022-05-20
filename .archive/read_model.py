# this script was originally meant to read multiple EBM models using the
# built-in web viewer (which we no longer use)

import sys
from time import sleep
import pickle
from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show


def show_model(path):
    with open(path, "rb") as f:
        ebm = pickle.load(f)

    ebm_global = ebm.explain_global()
    show(ebm_global)


for p in sys.argv[1:]:
    show_model(p)

while True:
    sleep(1)
