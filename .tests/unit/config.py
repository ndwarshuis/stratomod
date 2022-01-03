import unittest as ut
import sys
import inspect
from os.path import dirname, abspath

currentdir = dirname(abspath(inspect.getfile(inspect.currentframe())))
parentdir = dirname(dirname(currentdir))
sys.path.append(parentdir)

import workflow.scripts.common.config as conf


def assert_merge_noop(glb, lcl):
    assert conf.merge_dicts(glb, lcl) == glb


class ConfigTest(ut.TestCase):
    def test_noop(self):
        glb = {"a": 1, "b": {"c": 2}}
        assert_merge_noop(glb, {})

    def test_single_overwrite(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"a": 2}
        mrg = {"a": 2, "b": {"c": 2}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_overwrite_noop(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"a": 1}
        assert_merge_noop(glb, lcl)

    def test_single_addition(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"d": 2}
        mrg = {"a": 1, "b": {"c": 2}, "d": 2}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_addition_empty(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"d": {}}
        mrg = {"a": 1, "b": {"c": 2}, "d": {}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_addition_nonempty(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"d": {"e": 5}}
        mrg = {"a": 1, "b": {"c": 2}, "d": {"e": 5}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_addition_noop(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"d": None}
        assert_merge_noop(glb, lcl)

    def test_single_deletion(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": None}
        mrg = {"a": 1}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_nested_overwrite(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"c": 3}}
        mrg = {"a": 1, "b": {"c": 3}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_nested_overwrite_noop(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"c": 2}}
        assert_merge_noop(glb, lcl)

    def test_single_nested_addition(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"d": 4}}
        mrg = {"a": 1, "b": {"c": 2, "d": 4}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_nested_empty(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"d": {}}}
        mrg = {"a": 1, "b": {"c": 2, "d": {}}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_nested_nonempty(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"d": {"e": 5}}}
        mrg = {"a": 1, "b": {"c": 2, "d": {"e": 5}}}
        assert conf.merge_dicts(glb, lcl) == mrg

    def test_single_nested_addition_noop(self):
        glb = {"a": 1, "b": {"c": 2}}
        lcl = {"b": {"c": 2, "d": None}}
        assert_merge_noop(glb, lcl)

    def test_single_nested_deletion(self):
        glb = {"a": 1, "b": {"c": 2, "d": 3}}
        lcl = {"b": {"d": None}}
        mrg = {"a": 1, "b": {"c": 2}}
        assert conf.merge_dicts(glb, lcl) == mrg


if __name__ == "__main__":
    ut.main()
