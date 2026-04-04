"""
plank_cmb.py — Planck PR4 CMB likelihood for CAR framework
SCT Cosmology Series | DR JM NIPOK (2026)

NOTE ON FILENAME: This file is intentionally named plank_cmb.py (missing 'n')
for historical consistency with the repository. All imports use this spelling.

FIXES IN v2.0:
  - using_mock flag added for circular-data warning system
  - theta_star units clarified: this likelihood expects 100×θ* (e.g. 1.04105),
    NOT radians. The compute_chi2 comparison is in 100×θ* units throughout.
  - Planck 2018 central value confirmed: 100θ* = 1.04105 ± 0.00031
"""

import numpy as np
import os


class Planck_Likelihood:
    """
    Planck PR4 CMB angular scale likelihood.

    Operates on 100×θ* (the Planck-reported quantity).
    Planck 2018: 100θ* = 1.04105 ± 0.00031

    IMPORTANT: pass theta_star as 100×θ* (e.g. 1.04105), NOT in radians.
    The sct_core.py v2.0 CAR_predictions() dict returns theta_star=1.04105.
    """

    def __init__(self, data_dir=None):
        if data_dir is None:
            self.data_dir = os.path.join(
                os.path.dirname(__file__), '..', 'data')
        else:
            self.data_dir = data_dir

        self.ndata      = 86
        self.using_mock = False
        self.load_data()

    def load_data(self):
        planck_file = os.path.join(self.data_dir, 'planck_pr4_like.txt')

        if os.path.exists(planck_file):
            data = np.loadtxt(planck_file)
            self.theta_star_obs = float(data[0])   # in units of 100×θ*
            self.theta_star_err = float(data[1])
            self.using_mock = False
        else:
            import warnings
            warnings.warn(
                "Planck: using 2018 published value 100θ*=1.04105±0.00031. "
                "Run data/download_data.sh for full Planck PR4 likelihood.",
                UserWarning, stacklevel=2
            )
            self.using_mock = True
            # Planck 2018 Table 2 (Planck Collaboration 2020)
            self.theta_star_obs = 1.04105    # 100×θ*, NOT degrees, NOT radians
            self.theta_star_err = 0.00031

    def compute_chi2(self, theta_star_100):
        """
        Chi-squared for Planck angular scale.

        Parameters
        ----------
        theta_star_100 : float
            100×θ* (e.g. 1.04105 — the value Planck publishes)
            NOT radians, NOT degrees.
        """
        diff = (self.theta_star_obs - theta_star_100) / self.theta_star_err
        return diff ** 2

    def log_likelihood(self, theta_star_100):
        return -0.5 * self.compute_chi2(theta_star_100)
