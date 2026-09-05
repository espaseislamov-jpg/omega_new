import unittest
from pathlib import Path

import pandas as pd

from omega_core import io, matching, metrics, rt_profile


class LegacyC20MarkupTests(unittest.TestCase):
    def test_markup_uses_resolved_cluster_centres_not_placeholder_json_rts(self):
        project_dir = Path(__file__).resolve().parents[1]
        targets = io.load_reference_targets(project_dir / "reference_targets_reverted_c22fixed.json")
        retention_times = list(rt_profile.MANUAL_TABLE_RTS.values())
        peaks = pd.DataFrame({
            "peak_id": range(len(retention_times)),
            "apex_x": retention_times,
            "area": [1000.0] * len(retention_times),
            "percent_area": [100.0 / len(retention_times)] * len(retention_times),
            "prominence": [1000.0] * len(retention_times),
            "start_x": [value - 0.01 for value in retention_times],
            "end_x": [value + 0.01 for value in retention_times],
        })

        matched, shift = matching.match_targets_to_peaks(targets, peaks)
        corrected = matched.set_index("code")["corrected_target_rt"]

        self.assertAlmostEqual(corrected["C20:4N6"], 8.381 + shift)
        self.assertAlmostEqual(corrected["C20:5"], 8.410 + shift)
        self.assertAlmostEqual(corrected["C20:3N8"], 8.467 + shift)
        self.assertGreater(corrected["C20:3N8"], corrected["C20:5"])

    def test_display_markup_does_not_retune_legacy_judge(self):
        project_dir = Path(__file__).resolve().parents[1]
        targets = io.load_reference_targets(project_dir / "reference_targets_reverted_c22fixed.json")
        retention_times = list(rt_profile.MANUAL_TABLE_RTS.values())
        peaks = pd.DataFrame({
            "peak_id": range(len(retention_times)),
            "apex_x": retention_times,
            "area": [1000.0] * len(retention_times),
            "percent_area": [100.0 / len(retention_times)] * len(retention_times),
            "prominence": [1000.0] * len(retention_times),
            "start_x": [value - 0.01 for value in retention_times],
            "end_x": [value + 0.01 for value in retention_times],
        })
        matched, _ = matching.match_targets_to_peaks(targets, peaks)
        omega = metrics.compute_omega(matched)
        original_risk = metrics.assess_high_error_risk(matched, omega)

        changed_display = matched.copy()
        row_idx = changed_display.index[changed_display["code"] == "C20:3N8"][0]
        changed_display.at[row_idx, "corrected_target_rt"] += 1.0
        changed_risk = metrics.assess_high_error_risk(changed_display, omega)

        self.assertEqual(original_risk, changed_risk)


if __name__ == "__main__":
    unittest.main()
