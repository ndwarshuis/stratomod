from more_itertools import flatten


def lookup_config(config, *keys):
    k = keys[0]
    ks = keys[1:]
    return config[k] if len(ks) == 0 else lookup_config(config[k], *ks)


def lookup_global_chr_filter(config):
    return [
        *set(flatten(i["chr_filter"] for i in config["inputs"].values() if len(i) > 0))
    ]
