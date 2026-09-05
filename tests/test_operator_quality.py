from __future__ import annotations

import unittest

import pandas as pd
from matplotlib.figure import Figure

from New_idea import ChromatogramApp


class OperatorQualityTests(unittest.TestCase):
    def test_plot_view_can_be_restored_after_redraw(self):
        figure = Figure()
        axes = [figure.add_subplot(211), figure.add_subplot(212)]
        axes[0].set_xlim(8.31, 8.48)
        axes[0].set_ylim(-0.02, 1.15)
        axes[1].set_xlim(9.21, 9.32)
        axes[1].set_ylim(0.10, 0.92)
        saved = ChromatogramApp._capture_axes_views(axes)

        for axis in axes:
            axis.clear()
        ChromatogramApp._restore_axes_views(axes, saved)

        self.assertEqual(tuple(axes[0].get_xlim()), (8.31, 8.48))
        self.assertEqual(tuple(axes[0].get_ylim()), (-0.02, 1.15))
        self.assertEqual(tuple(axes[1].get_xlim()), (9.21, 9.32))
        self.assertEqual(tuple(axes[1].get_ylim()), (0.10, 0.92))

    def test_only_strong_geometry_is_green(self):
        good = ChromatogramApp._operator_quality({
            "geometry_score": 85.0,
            "high_error_risk": {"score": 0},
        })
        check = ChromatogramApp._operator_quality({
            "geometry_score": 84.9,
            "high_error_risk": {"score": 0},
        })
        stop = ChromatogramApp._operator_quality({
            "geometry_score": 54.9,
            "high_error_risk": {"score": 0},
        })

        self.assertEqual(good[1], "quality_good")
        self.assertEqual(check[1], "quality_check")
        self.assertEqual(stop[1], "quality_stop")

    def test_high_error_rule_still_has_priority(self):
        result = ChromatogramApp._operator_quality({
            "geometry_score": 95.0,
            "high_error_risk": {
                "score": 88,
                "peak_codes": ["C22:5", "C22:4"],
            },
        })
        self.assertEqual(result[1], "quality_check")

    def test_batch_export_places_requested_acids_after_quality(self):
        app = object.__new__(ChromatogramApp)
        app.loaded_batches = [{
            "sample_name": "O1_12345.D",
            "omega_report": 5.25,
            "confidence": {"geometry_score": 90.0, "high_error_risk": {"score": 0}},
            "matched_targets_df": pd.DataFrame([
                {"code": "C20:4N6", "percent_area": 8.1},
                {"code": "C18:2N6C", "percent_area": 12.2},
                {"code": "C18:3N6", "percent_area": 0.35},
            ]),
        }]

        row = app.build_batch_results_rows()[0]

        self.assertEqual(row[1:7], (
            "O1_12345.D",
            "5.2500",
            "ГОТОВО — правка не нужна",
            "8.1000",
            "12.2000",
            "0.3500",
        ))


if __name__ == "__main__":
    unittest.main()
