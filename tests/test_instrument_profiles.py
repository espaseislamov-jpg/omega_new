from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import omega_core
from omega_core import instrument_profiles
from omega_core import matching, metrics, rt_profile


class InstrumentProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.targets = omega_core.load_reference_targets(Path("reference_targets_reverted_c22fixed.json"))

    def test_legacy_profile_preserves_reference_targets(self):
        profile = instrument_profiles.legacy_profile(self.targets)
        applied = instrument_profiles.apply_profile_to_targets(self.targets, profile)
        self.assertTrue(applied["expected_rt"].equals(self.targets["expected_rt"]))
        self.assertTrue(applied["rt_reliable"].equals(self.targets["rt_reliable"]))
        self.assertFalse(applied["instrument_profile_custom_rt"].any())
        self.assertTrue((applied["result_multiplier"] == 1.0).all())
        self.assertTrue(profile["judge_calibrated"])
        self.assertTrue(applied["instrument_profile_judge_calibrated"].all())
        self.assertAlmostEqual(
            float(applied["instrument_profile_c22_height_priority_threshold"].iloc[0]),
            instrument_profiles.LEGACY_C22_HEIGHT_PRIORITY_THRESHOLD,
        )

    def test_profile_store_round_trip_and_active_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profiles.json"
            store = instrument_profiles.default_store(self.targets)
            custom = instrument_profiles.new_profile("GC-2", self.targets)
            custom["result_multiplier"] = 0.97
            store["profiles"].append(custom)
            store["active_profile_id"] = custom["id"]
            instrument_profiles.save_store(store, self.targets, path)
            loaded = instrument_profiles.load_store(self.targets, path)
        self.assertEqual(instrument_profiles.active_profile(loaded)["name"], "GC-2")
        self.assertAlmostEqual(instrument_profiles.active_profile(loaded)["result_multiplier"], 0.97)
        self.assertEqual(loaded["profiles"][0]["id"], instrument_profiles.LEGACY_PROFILE_ID)
        self.assertFalse(instrument_profiles.active_profile(loaded)["judge_calibrated"])
        self.assertEqual(
            instrument_profiles.judge_calibration_label(instrument_profiles.active_profile(loaded)),
            "Не настроен",
        )

    def test_export_import_creates_independent_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            store = instrument_profiles.default_store(self.targets)
            profile = instrument_profiles.new_profile("Прибор лаборатории", self.targets)
            profile["result_multiplier"] = 1.125
            profile["calculation_mode"] = instrument_profiles.DIRECT_COMPONENTS_CALCULATION_MODE
            profile["omega_component_codes"] = ["C20:3N8", "C22:6", "C22:4"]
            export_path = Path(directory) / "lab.omega-profile.json"
            instrument_profiles.export_profile(profile, export_path, self.targets["code"])
            raw = json.loads(export_path.read_text(encoding="utf-8"))
            imported = instrument_profiles.import_profile(export_path, store, self.targets)
        self.assertEqual(raw["schema_version"], 1)
        self.assertNotEqual(imported["id"], profile["id"])
        self.assertEqual(imported["name"], profile["name"])
        self.assertAlmostEqual(imported["result_multiplier"], 1.125)
        self.assertEqual(
            imported["calculation_mode"],
            instrument_profiles.DIRECT_COMPONENTS_CALCULATION_MODE,
        )
        self.assertEqual(imported["omega_component_codes"], ["C20:3N8", "C22:6", "C22:4"])
        self.assertEqual(imported["retention_times"], profile["retention_times"])
        self.assertFalse(imported["judge_calibrated"])

    def test_custom_profile_requires_strictly_ordered_times(self):
        profile = instrument_profiles.new_profile("Bad", self.targets)
        profile["retention_times"]["C20:5"] = profile["retention_times"]["C20:4N6"]
        with self.assertRaisesRegex(ValueError, "строго возрастать"):
            instrument_profiles.validate_profile(profile, self.targets["code"])

    def test_custom_profile_replaces_all_target_times(self):
        profile = instrument_profiles.new_profile("Shifted", self.targets)
        profile["retention_times"] = {
            code: value + 0.12 for code, value in profile["retention_times"].items()
        }
        applied = instrument_profiles.apply_profile_to_targets(self.targets, profile)
        expected = applied.set_index("code")["expected_rt"].to_dict()
        self.assertTrue(applied["rt_reliable"].all())
        self.assertAlmostEqual(expected["C20:5"], profile["retention_times"]["C20:5"])
        self.assertTrue(applied["instrument_profile_custom_rt"].all())

    def test_multiplier_changes_only_reported_omega_values(self):
        result = {
            "omega": {
                "omega3_trio": 4.0,
                "omega3_trio_strict": 3.8,
                "omega3_trio_corrected": 4.0,
                "total_area": 1234.0,
            },
            "omega_report": 4.0,
            "confidence": {"score": 91.0},
        }
        profile = {"id": "x", "name": "Column B", "result_multiplier": 1.1}
        scaled = instrument_profiles.apply_result_multiplier(result, profile)
        self.assertAlmostEqual(scaled["omega_report"], 4.4)
        self.assertAlmostEqual(scaled["omega"]["omega3_trio_strict"], 4.18)
        self.assertAlmostEqual(scaled["omega"]["omega3_trio_unscaled"], 4.0)
        self.assertAlmostEqual(scaled["omega"]["total_area"], 1234.0)
        self.assertEqual(scaled["confidence"], result["confidence"])
        self.assertTrue(np.isfinite(scaled["omega_report"]))

    def test_direct_component_profile_uses_selected_numerator(self):
        profile = instrument_profiles.new_profile("FID-2", self.targets)
        profile["calculation_mode"] = instrument_profiles.DIRECT_COMPONENTS_CALCULATION_MODE
        profile["omega_component_codes"] = ["C20:3N8", "C22:6", "C22:4"]
        profile = instrument_profiles.validate_profile(profile, self.targets["code"])
        applied = instrument_profiles.apply_profile_to_targets(self.targets, profile)

        areas = {code: 100.0 for code in applied["code"]}
        areas.update({"C20:3N8": 20.0, "C22:6": 30.0, "C22:4": 10.0})
        applied["area"] = applied["code"].map(areas).astype(float)
        applied["found_rt"] = applied["expected_rt"]
        applied["integration_start_x"] = applied["found_rt"] - 0.01
        applied["integration_end_x"] = applied["found_rt"] + 0.01
        applied["status"] = "matched_rt"

        result = metrics.compute_omega(applied)
        expected = 100.0 * 60.0 / applied["area"].sum()
        self.assertAlmostEqual(result["omega3_trio"], expected)
        self.assertAlmostEqual(result["omega3_trio_strict"], expected)
        self.assertTrue(result["profile_calculation_applied"])
        self.assertEqual(result["epa_overlap_credit_area"], 0.0)
        self.assertEqual(result["c22_overlap_credit_area"], 0.0)
        self.assertEqual(result["c22_overintegration_debit_points"], 0.0)
        risk = metrics.assess_high_error_risk(applied, result)
        self.assertGreaterEqual(risk["score"], 85)
        self.assertIn("instrument_profile_uncalibrated", risk["reason_codes"])

    def test_emergency_judge_uses_profile_components_and_relative_geometry(self):
        profile = instrument_profiles.new_profile("FID-2 emergency", self.targets)
        profile["calculation_mode"] = instrument_profiles.DIRECT_COMPONENTS_CALCULATION_MODE
        profile["omega_component_codes"] = ["C20:3N8", "C22:6", "C22:4"]
        profile["judge_emergency_enabled"] = True
        profile["judge_target_abs_error"] = 0.8
        profile = instrument_profiles.validate_profile(profile, self.targets["code"])
        applied = instrument_profiles.apply_profile_to_targets(self.targets, profile)

        applied["area"] = 100.0
        applied["found_rt"] = applied["expected_rt"]
        applied["corrected_target_rt"] = applied["expected_rt"]
        applied["integration_start_x"] = applied["found_rt"] - 0.015
        applied["integration_end_x"] = applied["found_rt"] + 0.015
        applied["matched_peak_id"] = np.arange(len(applied), dtype=float)
        applied["peak_height_smooth"] = 100.0
        applied.loc[applied["code"] == "C20:5", "peak_height_smooth"] = 50.0
        applied["status"] = "matched_rule"

        omega = metrics.compute_omega(applied)
        passed = metrics.assess_high_error_risk(applied, omega)
        self.assertEqual(passed["score"], 0)
        self.assertTrue(passed["instrument_profile_judge_emergency_enabled"])

        shifted = applied.copy()
        shifted.loc[shifted["code"] == "C20:3N8", "found_rt"] += 0.03
        warned = metrics.assess_high_error_risk(shifted, metrics.compute_omega(shifted))
        self.assertGreaterEqual(warned["score"], 85)
        self.assertIn("emergency_profile_geometry", warned["reason_codes"])
        self.assertIn("C20:3N8", warned["peak_codes"])

        missing = applied.copy()
        missing.loc[missing["code"] == "C20:3N8", "area"] = np.nan
        stopped = metrics.assess_high_error_risk(missing, metrics.compute_omega(missing))
        self.assertGreaterEqual(stopped["score"], 95)
        self.assertIn("C20:3N8", stopped["peak_codes"])

    def test_emergency_judge_metadata_survives_profile_round_trip(self):
        profile = instrument_profiles.new_profile("FID-2 emergency", self.targets)
        profile["judge_emergency_enabled"] = True
        profile["judge_target_abs_error"] = 0.8
        profile["judge_manual_samples"] = 19
        validated = instrument_profiles.validate_profile(profile, self.targets["code"])
        applied = instrument_profiles.apply_profile_to_targets(self.targets, validated)
        restored = instrument_profiles.profile_from_targets(applied)

        self.assertTrue(restored["judge_emergency_enabled"])
        self.assertAlmostEqual(restored["judge_target_abs_error"], 0.8)
        self.assertEqual(
            instrument_profiles.judge_calibration_label(restored),
            "Экстренный ±0.8, 19 проб",
        )

    def test_profile_calibration_follows_gradual_rt_drift(self):
        expected = np.asarray(list(rt_profile.MANUAL_TABLE_RTS.values()), dtype=float)
        observed = 1.004 * expected - 0.025
        observed = np.concatenate([observed, [6.81, 8.72, 10.12]])
        calibrated, slope, intercept = matching.estimate_profile_rt_calibration(expected, observed)
        self.assertAlmostEqual(slope, 1.004, places=5)
        self.assertAlmostEqual(intercept, -0.025, places=5)
        self.assertLess(float(np.max(np.abs(calibrated - (1.004 * expected - 0.025)))), 1e-8)

    def test_preview_windows_move_each_cluster_by_its_local_correction(self):
        profile = instrument_profiles.new_profile("Shifted", self.targets)
        increments = {
            "C16": 0.40,
            "C18": 0.45,
            "C20": 0.50,
            "C22": 0.55,
            "C24": 0.60,
        }
        profile["retention_times"] = {
            code: value + increments[code.split(":", 1)[0]]
            for code, value in profile["retention_times"].items()
        }
        applied = instrument_profiles.apply_profile_to_targets(self.targets, profile)
        windows = rt_profile.preview_windows(applied)
        self.assertAlmostEqual(windows[0][1], 6.40, places=6)
        self.assertAlmostEqual(windows[1][1], 7.85, places=6)
        self.assertAlmostEqual(windows[2][1], 8.80, places=6)
        self.assertAlmostEqual(windows[3][1], 9.65, places=6)


if __name__ == "__main__":
    unittest.main()
