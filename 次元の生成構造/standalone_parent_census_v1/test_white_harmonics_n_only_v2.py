#!/usr/bin/env python3
"""Nのみ白色零閉塞親生成器・独立粒子表の回帰試験。"""

from __future__ import annotations

import importlib.util
import hashlib
import json
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


generator = load_module(
    "make_parent_white_harmonics_n_only_v2",
    "make_parent_white_harmonics_n_only_v2.py",
)
reader = load_module(
    "particle_table_white_harmonics_n_only_v2",
    "particle_table_white_harmonics_n_only_v2.py",
)


class WhiteHarmonicParentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = generator.make_parent(5, seed=2, max_retries=99)

    def test_shape_and_hierarchical_zero_closure(self) -> None:
        result = self.result
        self.assertEqual(result.m, 5 * 4 // 2)
        self.assertEqual(result.relation_waves.shape, (10, 5))
        row_power = np.sum(np.abs(result.relation_waves) ** 2, axis=1)
        row_closure = np.sum(result.relation_waves**2, axis=1)
        self.assertLess(float(np.max(np.abs(row_closure) / row_power)), 1e-12)
        self.assertLess(
            float(abs(np.sum(result.relation_waves**2))),
            1e-12 * float(np.sum(row_power)),
        )
        self.assertLess(float(abs(result.parent_vector @ result.parent_vector)), 1e-12)

    def test_explicit_seed_is_bitwise_reproducible(self) -> None:
        repeated = generator.make_parent(5, seed=2, max_retries=0)
        self.assertTrue(np.array_equal(self.result.raw_white_noise, repeated.raw_white_noise))
        self.assertTrue(np.array_equal(self.result.parent_vector, repeated.parent_vector))
        self.assertTrue(np.array_equal(self.result.relation_waves, repeated.relation_waves))

    def test_explicit_seed_aborts_after_one_failed_attempt(self) -> None:
        failed = generator.SeedAttempt(
            attempt_number=0,
            seed=0,
            converged=False,
            iterations=10,
            residual=1.0,
            mu=0.0,
            sigma_max=0.0,
            parent_closure_abs=1.0,
            failure_reason="test_failure",
            parent_vector=None,
        )
        with mock.patch.object(generator, "_attempt_seed", return_value=(failed, None)) as call:
            with self.assertRaises(generator.ParentConstructionError) as caught:
                generator.make_parent(5, seed=123, max_retries=99)
        self.assertEqual(call.call_count, 1)
        self.assertEqual(len(caught.exception.attempts), 1)

    def test_manifest_has_no_harmonic_stage_parameter(self) -> None:
        manifest = generator.base_manifest(5, 2, 3)
        absent = manifest["explicit_absences"]
        self.assertIsNone(absent["harmonic_count_H"])
        self.assertFalse(absent["preassigned_harmonic_orders"])
        self.assertFalse(absent["preassigned_harmonic_amplitudes"])
        self.assertFalse(absent["preassigned_harmonic_phases"])
        self.assertEqual(manifest["derived_only_from_N"]["M"], 10)
        self.assertEqual(manifest["derived_only_from_N"]["sample_count_per_relation_wave"], 5)

    def test_reader_preserves_all_waves_phases_and_spectral_components(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent_dir = Path(temporary) / "parent"
            manifest = generator.base_manifest(5, 2, 3)
            manifest["seed_selection_protocol"] = (
                "integer seeds tested sequentially from 1; first converged seed accepted"
            )
            generator.write_success(parent_dir, self.result, manifest)
            census = reader.analyse(parent_dir)

        summary = census["summary"]
        self.assertEqual(summary["relation_wave_state_count"], 10)
        self.assertEqual(summary["phase_distinguished_wave_count"], 10)
        self.assertEqual(summary["total_DFT_components"], 10 * 5)
        self.assertEqual(len(census["waves"]), 10)
        self.assertEqual(len(census["harmonics"]), 50)
        self.assertTrue(summary["all_relation_waves_zero_closed"])
        self.assertLess(summary["DFT_reconstruction_max_error"], 1e-12)
        for row in census["waves"]:
            self.assertEqual(
                row["detected_non_dc_component_count"],
                row["base_wave_component_count"] + row["detected_harmonic_count"],
            )
            self.assertIn("mass_squared_type_gram_det", row)
            self.assertIn("charge_signed_from_rotation", row)
            self.assertIn(row["BFE_phase_network_readout"], {"B", "F", "E", "—"})
            self.assertIn("spin_type_readout", row)
            self.assertIn("correlation_lifetime_samples", row)
            self.assertIn("state_return_order", row)
            self.assertIn("quadratic_return_order", row)
            self.assertIn("cover_ratio", row)
            self.assertTrue(row["zero_closed"])

    def test_BFE_phase_network_has_no_external_resolution(self) -> None:
        same = np.array([1.0 + 0.0j, 2.0 + 0.0j, 3.0 + 0.0j])
        opposite = np.array([1.0 + 0.0j, -2.0 + 0.0j])
        intermediate = np.array([1.0 + 0.0j, 1.0j, -1.0 + 0.0j])
        self.assertEqual(reader.phase_network_readout(same, 1e-12)[0], "B")
        self.assertEqual(reader.phase_network_readout(opposite, 1e-12)[0], "F")
        self.assertEqual(reader.phase_network_readout(intermediate, 1e-12)[0], "E")

    def test_spin_cover_readout_on_pure_odd_and_even_waves(self) -> None:
        grid = np.arange(8)
        odd = np.exp(2j * np.pi * grid / 8)
        even = np.exp(4j * np.pi * grid / 8)
        odd_state, odd_quad = reader.cyclic_return_orders(odd)
        even_state, even_quad = reader.cyclic_return_orders(even)
        odd_half = reader.half_turn_cover(odd)
        even_half = reader.half_turn_cover(even)
        self.assertEqual((odd_state, odd_quad), (8, 4))
        self.assertEqual((even_state, even_quad), (4, 2))
        self.assertEqual(odd_half[-1], 2.0)
        self.assertEqual(even_half[-1], 1.0)
        self.assertIn(
            "2:1被覆確認",
            reader.spin_type("純奇数倍音（F型）", odd_half[-1]),
        )
        self.assertIn(
            "1:1被覆確認",
            reader.spin_type("純偶数倍音（B型）", even_half[-1]),
        )


class SavedArtifactTests(unittest.TestCase):
    def test_saved_N5_N40_are_complete_and_reproducible(self) -> None:
        generator_sha = hashlib.sha256(
            (HERE / "make_parent_white_harmonics_n_only_v2.py").read_bytes()
        ).hexdigest()
        reader_sha = hashlib.sha256(
            (HERE / "particle_table_white_harmonics_n_only_v2.py").read_bytes()
        ).hexdigest()
        combined = json.loads(
            (HERE / "particle_tables_white_harmonics_N5_N40_v3" / "summary.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(combined["reader_sha256"], reader_sha)
        combined_markdown = (
            HERE / "particle_tables_white_harmonics_N5_N40_v3" / "summary.md"
        ).read_text(encoding="utf-8")
        for heading in ("B/F/E", "電荷型 q", "質量型 μ", "スピン型"):
            self.assertIn(heading, combined_markdown)
        self.assertEqual(combined_markdown.count("| 5 | W"), 10)
        self.assertEqual(combined_markdown.count("| 40 | W"), 780)

        expected = [(5, 10, 2), (40, 780, 1)]
        for n, m, seed in expected:
            parent_dir = HERE / f"parent_white_harmonics_N{n}_v2"
            table_dir = HERE / "particle_tables_white_harmonics_N5_N40_v3" / f"N{n}"
            manifest = json.loads((parent_dir / "manifest.json").read_text(encoding="utf-8"))
            census = json.loads((table_dir / "census.json").read_text(encoding="utf-8"))
            waves = np.load(parent_dir / "relation_waves.npy")
            self.assertEqual(manifest["generator_sha256"], generator_sha)
            self.assertEqual(manifest["accepted_seed"], str(seed))
            self.assertEqual(waves.shape, (m, n))
            self.assertEqual(len(census["waves"]), m)
            self.assertEqual(len(census["harmonics"]), m * n)
            self.assertEqual(census["summary"]["phase_distinguished_wave_count"], m)
            self.assertTrue(census["summary"]["all_relation_waves_zero_closed"])
            self.assertTrue(all(row["zero_closed"] for row in census["waves"]))
            self.assertTrue(
                all(
                    np.isfinite(item["amplitude"])
                    and np.isfinite(item["phase_deg"])
                    for item in census["harmonics"]
                )
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                repeated = generator.make_parent(n, seed=seed, max_retries=3)
            self.assertTrue(np.array_equal(waves, repeated.relation_waves))


if __name__ == "__main__":
    unittest.main(verbosity=2)
