from functools import reduce
from typing import TypeVar, Optional, Callable


X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


def compose(*fs: Callable) -> Callable:
    return reduce(lambda f, g: lambda x: f(g(x)), fs, lambda x: x)


def maybe(default: Y, f: Callable[[X], Y], x: Optional[X]) -> Y:
    return default if x is None else f(x)


def flip(f: Callable[[X, Y], Z]) -> Callable[[Y, X], Z]:
    return lambda y, x: f(x, y)
