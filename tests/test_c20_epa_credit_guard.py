import unittest

import numpy as np
import pandas as pd

from omega_core import metrics


def _matched_c20(neighbor_width: float) -> pd.DataFrame:
    rows = [
        ("C20:4N6", 2300.0, 8.383, 0.052, "matched_c20_fit"),
        ("C20:5", 40.334076, 8.415, 0.02038, "matched_c20_fit"),
        ("C20:3N8", 495.606117, 8.471, neighbor_width, "matched_c20_fit"),
        ("C22:6", 600.562489, 9.240, 0.035, "matched_rt"),
        ("C22:5", 302.615334, 9.273, 0.030, "matched_rt"),
        ("C22:4", 900.0, 9.302, 0.040, "matched_rt"),
        ("C18:3N3", 300.0, 7.653, 0.035, "matched_rt"),
        ("other", 20714.430735, 10.0, 0.040, "matched_rt"),
    ]
    frame = pd.DataFrame(rows, columns=["code", "area", "found_rt", "width", "status"])
    frame["integration_start_x"] = frame["found_rt"] - frame["width"] / 2.0
    frame["integration_end_x"] = frame["found_rt"] + frame["width"] / 2.0
    frame["matched_peak_id"] = np.arange(len(frame), dtype=float)
    return frame


class C20EpaCreditGuardTests(unittest.TestCase):
    def test_compact_neighbor_keeps_severe_underfit_floor(self):
        result = metrics.compute_omega(_matched_c20(0.0673))

        self.assertAlmostEqual(result["epa_overlap_fraction"], 0.65, places=6)

    def test_broad_neighbor_does_not_receive_sixty_five_percent_floor(self):
        result = metrics.compute_omega(_matched_c20(0.0792))

        self.assertLess(result["epa_overlap_fraction"], 0.65)
        self.assertLess(
            result["epa_overlap_credit_area"],
            0.65 * result["epa_neighbor_area"],
        )


if __name__ == "__main__":
    unittest.main()
