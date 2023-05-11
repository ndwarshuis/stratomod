import scripts.python.common.config as cfg


# ensure all wildcards get parsed correctly
wildcard_constraints:
    **cfg._constraints,


# declare a function I use everywhere
def expand_refkey_from_refsetkey(path, wildcards):
    return expand(
        path,
        allow_missing=True,
        ref_key=config.refsetkey_to_refkey(wildcards.refset_key),
    )
