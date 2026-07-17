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


if __name__ == "__main__":
    unittest.main()
