import pickle
from dash import html
from sklearn import metrics
from sklearn.model_selection import train_test_split
from interpret import set_visualize_provider
from interpret.provider import InlineProvider
from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

# TODO some of these imports probably aren't actually needed


def write_model(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def read_model(path):
    with open(path, "rb") as f:
        return pickle.load(f)
