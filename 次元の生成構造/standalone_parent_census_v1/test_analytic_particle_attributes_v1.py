#!/usr/bin/env python3

import importlib.util
import math
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("analytic_particle_attributes_v1.py")
SPEC = importlib.util.spec_from_file_location("analytic_particle_attributes_v1", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class AnalyticParticleAttributesTest(unittest.TestCase):
    def test_allowed_modes_n5(self) -> None:
        self.assertEqual(MODULE.divisors(5), [1, 5])
        modes = MODULE.harmonic_modes(5)
        self.assertEqual(
            [(m.wavelength_multiple, m.harmonic_order, m.parity) for m in modes],
            [(1, 5, "奇数")],
        )

    def test_counts_n5(self) -> None:
        families, _ = MODULE.build_families(5)
        self.assertEqual(len(families), 3)
        self.assertEqual(sum(row.waveform_count for row in families), 8)
        self.assertEqual(MODULE.expected_waveform_count(5), 8)
        self.assertEqual(len(MODULE.exact_waveforms(5)), 8)

    def test_counts_n2_n3_n4(self) -> None:
        expected = {
            2: (3, 8, 0),
            3: (3, 8, 4),
            4: (6, 26, 0),
        }
        for n, (family_count, waveform_count, locked_count) in expected.items():
            with self.subTest(n=n):
                families, _ = MODULE.build_families(n)
                self.assertEqual(len(families), family_count)
                self.assertEqual(sum(row.waveform_count for row in families), waveform_count)
                self.assertEqual(
                    sum(
                        row.waveform_count
                        for row in families
                        if row.finite_order_lock and row.reflection_value > 0
                    ),
                    locked_count,
                )

    def test_counts_n40(self) -> None:
        families, _ = MODULE.build_families(40)
        self.assertEqual(MODULE.divisors(40), [1, 2, 4, 5, 8, 10, 20, 40])
        self.assertEqual(len(families), 42)
        self.assertEqual(sum(row.waveform_count for row in families), 5000)
        self.assertEqual(MODULE.expected_waveform_count(40), 5000)
        self.assertEqual(
            sum(
                row.waveform_count
                for row in families
                if row.finite_order_lock and row.reflection_value > 0
            ),
            112,
        )

    def test_reflection_and_spin(self) -> None:
        self.assertEqual(MODULE.reflection(1, 0), ("1/2", 0.5))
        self.assertEqual(MODULE.reflection(0, 2), ("0", 0.0))
        self.assertEqual(MODULE.reflection(1, 1), ("1/4", 0.25))
        self.assertEqual(MODULE.reflection(1, 2), ("1/6", 1 / 6))
        self.assertIn("2:1", MODULE.spin_readout(1, 0))
        self.assertIn("1:1", MODULE.spin_readout(0, 1))
        self.assertIn("混合", MODULE.spin_readout(1, 1))

    def test_address_charge(self) -> None:
        theta, finite, period, address, conjugate, charge, directions, _ = (
            MODULE.address_and_charge("1/2")
        )
        self.assertEqual(theta, "1/4")
        self.assertTrue(finite)
        self.assertEqual(period, 8)
        self.assertEqual(address, "+1/4")
        self.assertIn("3/4", conjugate)
        self.assertTrue(math.isclose(charge, 0.5))
        self.assertEqual(directions, 2)

        _, finite, _, address, _, charge, _, _ = MODULE.address_and_charge("1/6")
        self.assertFalse(finite)
        self.assertIn("非成立", address)
        self.assertIsNone(charge)


if __name__ == "__main__":
    unittest.main()
