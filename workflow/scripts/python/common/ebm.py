import pickle

# from dash import html  # type: ignore
# from sklearn import metrics  # type: ignore
# from sklearn.model_selection import train_test_split  # type: ignore
# from interpret import set_visualize_provider  # type: ignore
# from interpret.provider import InlineProvider  # type: ignore
from interpret.glassbox import ExplainableBoostingClassifier  # type: ignore

# from interpret import show  # type: ignore

# TODO some of these imports probably aren't actually needed


def write_model(path: str, obj: ExplainableBoostingClassifier) -> None:
    with open(path, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def read_model(path: str) -> ExplainableBoostingClassifier:
    with open(path, "rb") as f:
        return pickle.load(f)
