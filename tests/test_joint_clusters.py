from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omega_core.joint_clusters import (
    C22_CODES,
    build_joint_c22_candidate,
    judge_joint_c22_candidate,
    refine_c22_jointly,
)


class JointClusterTests(unittest.TestCase):
    def _fixture(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        x = np.linspace(9.18, 9.35, 1021)
        y = (
            900.0 * np.exp(-0.5 * ((x - 9.225) / 0.010) ** 2)
            + 460.0 * np.exp(-0.5 * ((x - 9.260) / 0.009) ** 2)
            + 280.0 * np.exp(-0.5 * ((x - 9.295) / 0.010) ** 2)
        )
        processed = pd.DataFrame({"x_corrected": x, "y_corrected": y, "y_smooth": y})
        rows = []
        for code, rt, start, end in [
            ("C22:6", 9.225, 9.190, 9.246),
            ("C22:5", 9.260, 9.246, 9.278),
            ("C22:4", 9.295, 9.278, 9.335),
        ]:
            mask = (x >= start) & (x <= end)
            rows.append({
                "code": code,
                "found_rt": rt,
                "area": float(np.trapezoid(y[mask], x[mask])),
                "integration_start_x": start,
                "integration_end_x": end,
                "percent_area": 1.0,
                "status": "synthetic",
            })
        return processed, pd.DataFrame(rows)

    def test_joint_partition_shares_valleys_and_conserves_cluster_area(self):
        processed, matched = self._fixture()
        old_total = float(matched["area"].sum())

        refined, decision = refine_c22_jointly(processed, matched)

        self.assertTrue(decision["accepted"])
        refined = refined.set_index("code")
        self.assertAlmostEqual(
            float(refined.at["C22:6", "integration_end_x"]),
            float(refined.at["C22:5", "integration_start_x"]),
        )
        self.assertAlmostEqual(
            float(refined.at["C22:5", "integration_end_x"]),
            float(refined.at["C22:4", "integration_start_x"]),
        )
        self.assertAlmostEqual(float(refined.loc[list(C22_CODES), "area"].sum()), old_total, places=9)

    def test_incomplete_cluster_is_rejected_without_partial_changes(self):
        processed, matched = self._fixture()
        incomplete = matched[matched["code"] != "C22:5"].copy()

        refined, decision = refine_c22_jointly(processed, incomplete)

        self.assertFalse(decision["accepted"])
        pd.testing.assert_frame_equal(refined, incomplete)

    def test_retry_is_limited_to_high_risk_and_small_result_shift(self):
        decision = {"accepted": True}
        candidate = {"omega_report": 4.20, "confidence": {"high_error_risk": {"score": 88}}}
        green = {"omega_report": 4.00, "confidence": {"high_error_risk": {"score": 0}}}
        red = {"omega_report": 4.00, "confidence": {"high_error_risk": {"score": 88}}}

        self.assertEqual(
            judge_joint_c22_candidate(green, candidate, decision)[1],
            "current_not_high_risk",
        )
        self.assertEqual(
            judge_joint_c22_candidate(red, candidate, decision),
            (True, "accepted"),
        )
        too_far = {"omega_report": 4.23, "confidence": {"high_error_risk": {"score": 88}}}
        self.assertEqual(
            judge_joint_c22_candidate(red, too_far, decision)[1],
            "omega_shift_too_large",
        )

    def test_custom_profile_is_not_repartitioned_without_validation(self):
        processed, matched = self._fixture()
        matched["instrument_profile_custom_rt"] = True

        candidate, decision = build_joint_c22_candidate({
            "processed_df": processed,
            "matched_targets_df": matched,
        })

        self.assertIsNone(candidate)
        self.assertEqual(decision["reason"], "custom_profile_not_validated")


if __name__ == "__main__":
    unittest.main()
