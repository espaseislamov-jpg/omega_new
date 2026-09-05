from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from omega_core import fit_recovery


class FitRecoverySafetyTests(unittest.TestCase):
    def test_invalid_start_is_rejected_without_calling_scipy(self):
        valid = fit_recovery._least_squares_inputs_are_valid(
            [np.nan, 1.0],
            [0.0, 0.0],
            [2.0, 2.0],
        )
        self.assertFalse(valid)

    def test_scipy_bounds_error_falls_back_instead_of_escaping(self):
        x = np.linspace(8.35, 8.50, 120)
        y = np.exp(-0.5 * ((x - 8.41) / 0.008) ** 2) * 100.0
        frame = pd.DataFrame({"x_corrected": x, "y_smooth": y})

        with patch.object(
            fit_recovery,
            "least_squares",
            side_effect=ValueError("Initial guess is outside of provided bounds"),
        ):
            result = fit_recovery._fit_cluster_components_gaussian(
                frame,
                initial_centers=[8.38, 8.41, 8.467],
                window_left=8.35,
                window_right=8.50,
                center_tolerances=[0.01, 0.01, 0.015],
            )

        self.assertIsNone(result)

    def test_missing_epa_fit_preserves_already_integrated_flanks(self):
        x = np.linspace(8.32, 8.50, 361)
        frame = pd.DataFrame({
            "x_corrected": x,
            "y_corrected": np.zeros_like(x),
            "y_smooth": np.zeros_like(x),
        })
        matched = pd.DataFrame([
            {
                "code": "C20:4N6", "found_rt": 8.357, "area": 3000.0,
                "integration_start_x": 8.330, "integration_end_x": 8.379,
                "matched_peak_id": 10.0, "status": "matched_c20_rule",
            },
            {
                "code": "C20:5", "found_rt": np.nan, "area": np.nan,
                "integration_start_x": 8.410, "integration_end_x": 8.466,
                "matched_peak_id": np.nan, "status": "not_found",
            },
            {
                "code": "C20:3N8", "found_rt": 8.443, "area": 700.0,
                "integration_start_x": 8.410, "integration_end_x": 8.466,
                "matched_peak_id": 11.0, "status": "matched_c20_rule",
            },
        ])
        peaks = pd.DataFrame({"peak_id": [10, 11], "apex_x": [8.357, 8.443]})
        fitted = [
            {"center": 8.357, "area": 1800.0, "amplitude": 1800.0, "fwhm": 0.014, "eta": 0.3},
            {"center": 8.389, "area": 14.0, "amplitude": 14.0, "fwhm": 0.014, "eta": 0.3},
            {"center": 8.443, "area": 320.0, "amplitude": 320.0, "fwhm": 0.014, "eta": 0.3},
        ]

        with patch.object(fit_recovery, "_fit_cluster_components", return_value=(fitted, {"r2": 0.995})):
            result = fit_recovery.recover_underintegrated_c20_components_with_fit(
                frame, peaks, matched
            ).set_index("code")

        self.assertEqual(result.at["C20:4N6", "area"], 3000.0)
        self.assertEqual(result.at["C20:3N8", "area"], 700.0)
        self.assertGreater(result.at["C20:5", "area"], 5.0)
        self.assertEqual(result.at["C20:5", "status"], "recovered_c20_missing_epa_component")
        self.assertGreaterEqual(result.at["C20:5", "integration_start_x"], 8.379)
        self.assertLessEqual(result.at["C20:5", "integration_end_x"], 8.410)

    def test_missing_epa_is_not_synthesized_without_both_flanks(self):
        frame = pd.DataFrame({
            "x_corrected": np.linspace(8.32, 8.50, 361),
            "y_corrected": np.zeros(361),
            "y_smooth": np.zeros(361),
        })
        matched = pd.DataFrame([
            {"code": "C20:4N6", "found_rt": 8.357, "area": 3000.0, "matched_peak_id": 10.0},
            {"code": "C20:5", "found_rt": np.nan, "area": np.nan, "matched_peak_id": np.nan},
            {"code": "C20:3N8", "found_rt": np.nan, "area": np.nan, "matched_peak_id": np.nan},
        ])
        peaks = pd.DataFrame({"peak_id": [10], "apex_x": [8.357]})

        with patch.object(fit_recovery, "_fit_cluster_components") as mocked_fit:
            result = fit_recovery.recover_underintegrated_c20_components_with_fit(
                frame, peaks, matched
            )

        self.assertFalse(mocked_fit.called)
        self.assertTrue(np.isnan(result.loc[result["code"] == "C20:5", "area"].iloc[0]))


if __name__ == "__main__":
    unittest.main()
