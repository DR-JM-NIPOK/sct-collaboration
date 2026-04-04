"""
des_y6_lensing.py — DES-Y6 weak lensing likelihood for CAR framework
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - Mock S8_obs changed from 0.783 (SCT prediction) to 0.780 ± 0.012
    which is the actual published DES-Y6 measurement
  - using_mock flag added for circular-data warning
  - Clear documentation of real data download path
"""

import numpy as np
import os


class DES_Y6_Likelihood:
    """
    DES-Y6 3×2pt weak lensing likelihood.

    Published value: S8 = 0.780 ± 0.012  (DES Collaboration 2026)

    IMPORTANT: The previous version used S8_obs = 0.783, which is SCT's
    own prediction. That made chi-squared ≈ 0 by construction and produced
    circular Bayesian evidence. This version uses the actual DES-Y6 value.
    """

    def __init__(self, data_dir=None):
        if data_dir is None:
            self.data_dir = os.path.join(
                os.path.dirname(__file__), '..', 'data')
        else:
            self.data_dir = data_dir

        self.ndata      = 1800
        self.using_mock = False
        self.load_data()

    def load_data(self):
        data_file = os.path.join(self.data_dir, 'des_y6_3x2pt.txt')

        if os.path.exists(data_file):
            data         = np.loadtxt(data_file)
            self.S8_obs  = float(data[0])
            self.S8_err  = float(data[1])
            self.using_mock = False
        else:
            import warnings
            warnings.warn(
                "DES-Y6: using published summary statistic S8=0.780±0.012 "
                "(DES Collaboration 2026). Run data/download_data.sh for "
                "full data vector.",
                UserWarning, stacklevel=2
            )
            self.using_mock = True
            # ACTUAL published DES-Y6 measurement — NOT SCT's prediction
            self.S8_obs = 0.780
            self.S8_err = 0.012

    def compute_chi2(self, S8):
        return ((self.S8_obs - S8) / self.S8_err) ** 2

    def log_likelihood(self, S8):
        return -0.5 * self.compute_chi2(S8)
