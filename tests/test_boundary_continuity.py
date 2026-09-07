import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

import omega_core
from omega_core import signal, instrument_profiles
from omega_core.boundary_audit import snapshot, record_changes
from New_idea import ChromatogramApp


class BoundaryContinuityTests(unittest.TestCase):
    def test_bright_narrow_target_is_recovered_without_optional_backend(self):
        x = np.linspace(7.60, 7.70, 1001)
        y = 100000 * np.exp(-0.5 * ((x - 7.650) / 0.002) ** 2)
        frame = pd.DataFrame({"x_corrected": x, "y_corrected": y,
                              "y_smooth": y, "dy": np.gradient(y, x)})
        with patch.object(signal, "ENABLE_PYOPENMS_PEAK_ASSIST", False):
            result = signal.detect_peak_candidates(frame)
        self.assertFalse(result.empty)
        self.assertLess(abs(result.iloc[0]["apex_x"] - 7.650), 0.001)
        self.assertGreater(result.iloc[0]["area"], 400)

    def test_supplement_does_not_replace_primary_peak_with_earlier_rt(self):
        primary = pd.DataFrame([{"apex_x": 1.0, "area": 100.0, "prominence": 100.0}])
        result = signal._merge_peak_records(primary, [
            {"apex_x": 0.997, "area": 2.0, "prominence": 2.0},
            {"apex_x": 2.0, "area": 50.0, "prominence": 50.0},
        ])
        self.assertEqual(result["apex_x"].tolist(), [1.0, 2.0])
        self.assertEqual(result["area"].tolist(), [100.0, 50.0])

    def test_empty_primary_detector_still_runs_recovery(self):
        x = np.linspace(0, 1, 1001)
        y = 1000 * np.exp(-0.5 * ((x - 0.5) / 0.02) ** 2)
        frame = pd.DataFrame({"x_corrected": x, "y_corrected": y, "y_smooth": y})
        recovered = pd.DataFrame({"peak_id": [1], "apex_x": [0.5]})
        with patch.object(signal, "find_peaks", return_value=(np.array([], dtype=int), {})), \
             patch.object(signal, "augment_targeted_cluster_peaks", return_value=recovered) as rescue, \
             patch.object(signal, "detect_peaks_with_pyopenms", return_value=pd.DataFrame()) as assist:
            result = signal.detect_peak_candidates(frame)
        rescue.assert_called_once()
        assist.assert_called_once()
        self.assertEqual(result.iloc[0]["apex_x"], 0.5)

    def test_area_rejection_still_runs_recovery(self):
        x = np.linspace(0, 1, 1001)
        y = 1000 * np.exp(-0.5 * ((x - 0.5) / 0.02) ** 2)
        frame = pd.DataFrame({"x_corrected": x, "y_corrected": y, "y_smooth": y})
        with patch.object(signal, "CHEMSTATION_INITIAL_AREA_REJECT", 1e20), \
             patch.object(signal, "augment_targeted_cluster_peaks", return_value=pd.DataFrame()) as rescue, \
             patch.object(signal, "detect_peaks_with_pyopenms", return_value=pd.DataFrame()):
            signal.detect_peak_candidates(frame)
        rescue.assert_called_once()

    def test_exact_profile_copy_preserves_engine_and_judge(self):
        targets = omega_core.load_reference_targets()
        source = instrument_profiles.legacy_profile(targets)
        copied = instrument_profiles.copy_profile(source, "Same instrument")
        self.assertNotEqual(source["id"], copied["id"])
        for key in source:
            if key not in {"id", "name"}:
                self.assertEqual(source[key], copied[key], key)
        store = instrument_profiles.normalize_store({"profiles": [copied]}, targets)
        persisted = next(p for p in store["profiles"] if p["id"] == copied["id"])
        self.assertTrue(persisted["judge_calibrated"])
        self.assertFalse(persisted["custom_rt"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.json"
            instrument_profiles.export_profile(copied, path, targets["code"])
            imported = instrument_profiles.import_profile(path, store, targets)
        self.assertFalse(imported["custom_rt"])
        self.assertTrue(imported["judge_calibrated"])
        code = next(iter(copied["retention_times"]))
        copied["retention_times"][code] += 1
        self.assertNotEqual(copied["retention_times"][code], source["retention_times"][code])

    def test_visual_bounds_are_numerical_bounds(self):
        row = {"integration_start_x": 1.0, "integration_end_x": 2.0, "found_rt": 1.5}
        self.assertEqual(ChromatogramApp._visual_peak_footprint_bounds(
            None, row, np.linspace(0, 3, 100), np.zeros(100)), (1.0, 2.0))

    def test_audit_records_removed_peak_without_changing_values(self):
        before = pd.DataFrame([{"code": "A", "area": 10.0, "found_rt": 1.0}])
        after = before.copy()
        after.loc[0, "area"] = np.nan
        expected = after.copy()
        history = []
        record_changes(snapshot(before), after, "test_guard", history)
        pd.testing.assert_frame_equal(after, expected)
        self.assertEqual(history[0]["before"]["area"], 10.0)
        self.assertIsNone(history[0]["after"]["area"])
        self.assertEqual(history[0]["stage"], "test_guard")
