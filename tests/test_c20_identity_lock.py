from __future__ import annotations

import unittest

import pandas as pd

from omega_core import clusters


def _matched(epa_rt: float, status: str = "matched_c20_rule") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "code": ["C20:4N6", "C20:5", "C20:3N8"],
            "found_rt": [8.375, epa_rt, 8.463],
            "area": [200.0, 40.0, 180.0],
            "matched_peak_id": [10, 11, 12],
            "match_score": [0.0, 0.0, 0.0],
            "status": ["matched", status, "matched"],
        }
    )


class C20IdentityLockTests(unittest.TestCase):
    def setUp(self):
        self.peaks = pd.DataFrame(
            {
                "peak_id": [10, 11, 13, 12],
                "apex_x": [8.375, 8.402, 8.418, 8.463],
            }
        )

    def test_fit_may_refine_inside_assigned_peak_cell(self):
        result = clusters._c20_identity_lock_after_fit(
            _matched(8.402),
            _matched(8.406, status="matched_c20_fit"),
            self.peaks,
        )

        self.assertAlmostEqual(result.loc[result.code == "C20:5", "found_rt"].iloc[0], 8.406)
        self.assertEqual(
            result.attrs[clusters.C20_ASSIGNMENT_TRACE_ATTR][-1]["decision"],
            "accepted",
        )

    def test_fit_cannot_cross_to_neighbor_peak(self):
        fitted = _matched(8.416, status="matched_c20_fit")
        fitted.loc[fitted.code == "C20:5", "area"] = 77.0
        result = clusters._c20_identity_lock_after_fit(
            _matched(8.402),
            fitted,
            self.peaks,
        )
        result = clusters._apply_c20_display_identity(result)
        row = result[result.code == "C20:5"].iloc[0]

        self.assertAlmostEqual(row.found_rt, 8.402)
        self.assertAlmostEqual(row.area, 77.0)
        self.assertIn("identity_center_locked", row.status)
        self.assertEqual(
            result.attrs[clusters.C20_ASSIGNMENT_TRACE_ATTR][-1]["decision"],
            "rejected_neighbor_jump",
        )

    def test_fit_hitting_center_tolerance_is_rejected_without_detected_neighbor(self):
        peaks = self.peaks[self.peaks.peak_id != 13].copy()
        result = clusters._c20_identity_lock_after_fit(
            _matched(8.402),
            _matched(8.412, status="matched_c20_fit"),
            peaks,
        )

        result = clusters._apply_c20_display_identity(result)
        row = result[result.code == "C20:5"].iloc[0]
        self.assertAlmostEqual(row.found_rt, 8.402)
        self.assertIn("identity_center_locked", row.status)


if __name__ == "__main__":
    unittest.main()
