import json


def merge_dicts(d1, d2):
    # NOTE: this will only recursively merge nested dicts (not lists, which
    # could also be present in a yaml config)
    # def use_other_maybe(k, v):
    def pick_dict(d1, d2, k):
        if k in d2:
            _v = d2[k]
            if isinstance(_v, dict):
                return merge_dicts(d1[k], _v) if k in d1 else _v
            else:
                return _v
        else:
            return d1[k]

    if d1 is None:
        return d2
    if d2 is None:
        return {}
    return {
        k: pick_dict(d1, d2, k)
        for k in set(d1.keys()) | set(d2.keys())
        if not (k in d2 and d2[k] is None)
    }


def lookup_run_config(config, run_key):
    r = config["ebm_runs"][run_key]
    g = config["global"]
    return merge_dicts(g, r)


def lookup_run_json(config, run_key):
    return json.dumps(lookup_run_config(config, run_key))


def lookup_config(config, *keys):
    k = keys[0]
    ks = keys[1:]
    return config[k] if len(ks) == 0 else lookup_config(config[k], *ks)
