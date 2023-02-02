from functools import reduce
from typing import TypeVar, Optional, Callable, Any


X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


def compose1(f1: Callable[[Y], Z], f2: Callable[[X], Y]) -> Callable[[X], Z]:
    return lambda x: f1(f2(x))


def compose(*fs: Callable[[Any], Any]) -> Callable[[Any], Any]:
    return reduce(compose1, fs, lambda x: x)


def maybe(default: Y, f: Callable[[X], Y], x: Optional[X]) -> Y:
    return default if x is None else f(x)


def flip(f: Callable[[X, Y], Z]) -> Callable[[Y, X], Z]:
    return lambda y, x: f(x, y)
