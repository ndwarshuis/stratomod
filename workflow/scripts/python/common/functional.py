from functools import reduce


def compose(*fs):
    return reduce(lambda f, g: lambda x: f(g(x)), fs, lambda x: x)
