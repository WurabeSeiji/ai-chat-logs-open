#!/usr/bin/env python3
"""make_parent_n_only_v1 の契約テスト。"""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("make_parent_n_only_tested", HERE / "make_parent_n_only_v1.py")
mp = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mp
SPEC.loader.exec_module(mp)


class MakeParentNOnlyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.first = mp.make_parent(5, seed=11, max_retries=99)
        cls.second = mp.make_parent(5, seed=11, max_retries=0)

    def test_shape_and_derived_relation_count(self) -> None:
        self.assertEqual(self.first.m, 10)
        self.assertEqual(self.first.parent_modes.shape, (10, 5))
        self.assertEqual(self.first.white_input.shape, (10, 5))

    def test_explicit_seed_reproducibility(self) -> None:
        self.assertTrue(np.array_equal(self.first.white_input, self.second.white_input))
        self.assertTrue(np.array_equal(self.first.parent_modes, self.second.parent_modes))
        self.assertEqual(len(self.first.attempts), 1)
        self.assertEqual(len(self.second.attempts), 1)

    def test_closure_and_white_weight_inheritance(self) -> None:
        closures = np.abs(np.sum(self.first.parent_modes**2, axis=0))
        self.assertLess(float(np.max(closures)), 1e-12)
        self.assertLess(float(abs(np.sum(self.first.parent_modes**2))), 1e-12)
        expected = np.linalg.norm(self.first.white_input, axis=0)
        expected /= np.linalg.norm(expected)
        self.assertTrue(np.array_equal(expected, self.first.mode_weights))
        self.assertAlmostEqual(float(np.linalg.norm(self.first.parent_modes)), 1.0, places=14)

    def test_explicit_unstable_seed_aborts_once(self) -> None:
        with self.assertRaises(mp.ParentConstructionError) as caught:
            mp.make_parent(5, seed=0, max_retries=99)
        self.assertEqual(len(caught.exception.attempts), 1)
        self.assertEqual(caught.exception.attempts[0].seed, 0)

    def test_internal_seed_retry_path(self) -> None:
        original = mp.secrets.randbits
        seeds = iter([0, 1, 11])
        mp.secrets.randbits = lambda bits: next(seeds)
        try:
            result = mp.make_parent(5, seed=None, max_retries=3)
        finally:
            mp.secrets.randbits = original
        self.assertEqual([attempt.seed for attempt in result.attempts], [0, 1, 11])
        self.assertEqual(result.seed, 11)
        self.assertFalse(result.seed_was_explicit)


if __name__ == "__main__":
    unittest.main(verbosity=2)
