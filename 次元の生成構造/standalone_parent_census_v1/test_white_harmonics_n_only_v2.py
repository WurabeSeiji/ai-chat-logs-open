#!/usr/bin/env python3
"""Nのみ白色零閉塞親・make_parent直後分類の回帰試験。"""

from __future__ import annotations

import importlib.util
import csv
import hashlib
import json
import sys
import tempfile
import unittest
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

    def test_reader_keeps_zero_closure_separate_from_stable_classification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent_dir = Path(temporary) / "parent"
            manifest = generator.base_manifest(5, 2, 3)
            manifest["seed_selection_protocol"] = (
                "integer seeds tested sequentially from 1; first converged seed accepted"
            )
            generator.write_success(parent_dir, self.result, manifest)
            census = reader.analyse(parent_dir)

        summary = census["summary"]
        self.assertEqual(summary["generated_relation_wave_count"], 10)
        self.assertEqual(summary["stable_classified_wave_count"], 0)
        self.assertEqual(summary["unclassified_wave_count"], 10)
        self.assertEqual(len(census["waves"]), 10)
        self.assertEqual(len(census["spectrum_by_cyclic_order"]), 5)
        self.assertTrue(summary["all_relation_waves_zero_closed"])
        self.assertLess(summary["DFT_reconstruction_max_error"], 1e-12)
        for row in census["waves"]:
            self.assertEqual(row["classification_status"], "分類外")
            self.assertIsNone(row["family_id"])
            self.assertGreater(row["outside_stationary_order_count"], 0)
            self.assertFalse(row["phase_is_classification_key"])
            self.assertFalse(row["phase_rounding_to_0_or_180"])
            self.assertTrue(row["zero_closed"])

    def test_n5_classification_uses_wavelength_not_phase(self) -> None:
        _, catalog = reader.load_catalog(5)
        first = np.zeros(5, dtype=complex)
        first[0] = 2.0 * np.exp(0.37j)  # lambda0のN点節上表現
        first[1] = 3.0 * np.exp(1.23j)  # 5lambda0
        second = first.copy()
        second[0] = 2.0 * np.exp(-2.41j)
        second[1] = 3.0 * np.exp(2.77j)

        one, _ = reader.classify_spectrum(first, 5, 10, catalog)
        two, _ = reader.classify_spectrum(second, 5, 10, catalog)
        self.assertEqual(one["family_id"], "N5-F003")
        self.assertEqual(two["family_id"], "N5-F003")
        self.assertEqual(one["present_harmonic_orders"], [5])
        self.assertFalse(one["phase_is_classification_key"])
        self.assertFalse(one["phase_rounding_to_0_or_180"])

    def test_n5_outside_wavelength_is_not_forced_into_family(self) -> None:
        _, catalog = reader.load_catalog(5)
        spectrum = np.zeros(5, dtype=complex)
        spectrum[0] = 1.0
        spectrum[1] = np.exp(0.4j)
        spectrum[2] = 0.5 * np.exp(-1.1j)  # 5/2 lambda0: 分類表外
        result, _ = reader.classify_spectrum(spectrum, 5, 10, catalog)
        self.assertEqual(result["classification_status"], "分類外")
        self.assertIsNone(result["family_id"])
        self.assertEqual(result["outside_stationary_orders"], [2])

    def test_n40_odd_harmonic_family_with_arbitrary_phases(self) -> None:
        _, catalog = reader.load_catalog(40)
        spectrum = np.zeros(40, dtype=complex)
        spectrum[4] = 2.0 * np.exp(0.83j)   # 10lambda0: 基底
        spectrum[20] = 0.7 * np.exp(-2.2j)  # 2lambda0: 5倍音
        result, _ = reader.classify_spectrum(spectrum, 40, 780, catalog)
        self.assertEqual(result["family_id"], "N40-F016")
        self.assertEqual(result["base_wavelength_multiple"], 10)
        self.assertEqual(result["present_harmonic_orders"], [5])
        self.assertAlmostEqual(result["actual_reflection_from_power"], 0.5)


class SavedArtifactTests(unittest.TestCase):
    def test_saved_N5_N40_immediate_classification_is_current(self) -> None:
        generator_sha = hashlib.sha256(
            (HERE / "make_parent_white_harmonics_n_only_v2.py").read_bytes()
        ).hexdigest()
        reader_sha = hashlib.sha256(
            (HERE / "particle_table_white_harmonics_n_only_v2.py").read_bytes()
        ).hexdigest()
        output = HERE / "stable_wave_classification_after_parent_N5_N40_N300_v1"
        combined = json.loads(
            (output / "summary.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(combined["reader_sha256"], reader_sha)
        self.assertEqual(combined["schema"], reader.OUTPUT_SCHEMA)
        self.assertEqual(
            [item["summary"]["N"] for item in combined["results"]], [5, 40, 300]
        )
        self.assertEqual(combined["results"][2]["summary"]["M"], 44850)

        expected = [(5, 10, 2), (40, 780, 1)]
        for n, m, seed in expected:
            parent_dir = HERE / f"parent_white_harmonics_N{n}_v2"
            table_dir = output / f"N{n}"
            manifest = json.loads((parent_dir / "manifest.json").read_text(encoding="utf-8"))
            census = json.loads((table_dir / "census.json").read_text(encoding="utf-8"))
            waves = np.load(parent_dir / "relation_waves.npy")
            self.assertEqual(manifest["generator_sha256"], generator_sha)
            self.assertEqual(manifest["accepted_seed"], str(seed))
            self.assertEqual(waves.shape, (m, n))
            self.assertEqual(census["wave_table"]["row_count"], m)
            with (table_dir / "make_parent直後の波分類.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                wave_rows = list(csv.DictReader(handle))
            self.assertEqual(len(wave_rows), m)
            self.assertEqual(len(census["spectrum_by_cyclic_order"]), n)
            self.assertEqual(census["summary"]["stable_classified_wave_count"], 0)
            self.assertEqual(census["summary"]["unclassified_wave_count"], m)
            self.assertTrue(census["summary"]["all_relation_waves_zero_closed"])
            self.assertTrue(all(row["安定波分類"] == "分類外" for row in wave_rows))
            self.assertTrue(
                all(
                    row["零閉塞相対残差_監査のみ"] not in {"", "—"}
                    for row in wave_rows
                )
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
