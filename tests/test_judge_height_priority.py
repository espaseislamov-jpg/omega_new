from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omega_core import metrics


class JudgeHeightPriorityTests(unittest.TestCase):
    def test_smoothed_peak_height_uses_only_accepted_integration_bounds(self):
        processed = pd.DataFrame({
            "x_corrected": np.arange(7, dtype=float),
            "y_smooth": [0.0, 8.0, 2.0, 50.0, 1.0, 4.0, 0.0],
        })
        matched = pd.DataFrame([
            {"code": "C22:5", "integration_start_x": 0.0, "integration_end_x": 2.0},
            {"code": "C22:4", "integration_start_x": 4.0, "integration_end_x": 6.0},
        ])

        annotated = metrics.annotate_peak_heights(processed, matched).set_index("code")

        self.assertAlmostEqual(float(annotated.at["C22:5", "peak_height_smooth"]), 8.0)
        self.assertAlmostEqual(float(annotated.at["C22:4", "peak_height_smooth"]), 4.0)

    def test_height_ratio_prioritizes_existing_warning_without_creating_one(self):
        low_risk = metrics.classify_high_error_risk({
            "instrument_profile_judge_calibrated": True,
            "C22_height_ratio": 2.0,
            "C22_height_priority_threshold": 1.3635,
        })
        warned = metrics.classify_high_error_risk({
            "instrument_profile_judge_calibrated": True,
            "C22_height_ratio": 2.0,
            "C22_height_priority_threshold": 1.3635,
            "omega_c22_overlap_legacy_fraction": 0.1,
            "C22_4_left_width": 0.02,
        })

        self.assertEqual(low_risk["score"], 0)
        self.assertEqual(low_risk["review_priority"], "normal")
        self.assertGreaterEqual(warned["score"], 85)
        self.assertEqual(warned["review_priority"], "c22_first")
        self.assertIn("c22_height_priority", warned["reason_codes"])

    def test_uncalibrated_custom_profile_cannot_be_green(self):
        risk = metrics.classify_high_error_risk({
            "instrument_profile_custom_rt": True,
            "instrument_profile_judge_calibrated": False,
        })

        self.assertGreaterEqual(risk["score"], 85)
        self.assertIn("instrument_profile_uncalibrated", risk["reason_codes"])
        self.assertFalse(risk["instrument_profile_judge_calibrated"])

    def test_emergency_custom_profile_can_pass_broad_geometry_gate(self):
        risk = metrics.classify_high_error_risk({
            "instrument_profile_custom_rt": True,
            "instrument_profile_judge_calibrated": False,
            "instrument_profile_judge_emergency_enabled": True,
            "instrument_profile_judge_target_abs_error": 0.8,
        })

        self.assertEqual(risk["score"], 0)
        self.assertNotIn("instrument_profile_uncalibrated", risk["reason_codes"])
        self.assertTrue(risk["instrument_profile_judge_emergency_enabled"])
        self.assertAlmostEqual(risk["instrument_profile_judge_target_abs_error"], 0.8)

    def test_emergency_geometry_warning_names_the_suspicious_peak(self):
        risk = metrics.classify_high_error_risk({
            "instrument_profile_custom_rt": True,
            "instrument_profile_judge_calibrated": False,
            "instrument_profile_judge_emergency_enabled": True,
            "instrument_profile_judge_target_abs_error": 0.8,
            "emergency_geometry_codes": ["C20:3N8"],
            "emergency_geometry_reasons": ["Вершина C20:3N8 смещена."],
        })

        self.assertGreaterEqual(risk["score"], 85)
        self.assertIn("emergency_profile_geometry", risk["reason_codes"])
        self.assertEqual(risk["peak_codes"], ["C20:3N8"])


if __name__ == "__main__":
    unittest.main()
