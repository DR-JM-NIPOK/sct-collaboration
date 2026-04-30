"""
combined_likelihood.py — Combined likelihood for all CAR datasets
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - Import corrected from planck_cmb → plank_cmb (matches actual filename)
  - Key names updated to match corrected sct_core.py v2.0 return dict
  - theta_star now passed correctly as 100×θ* value (1.04105), not radians
  - Added explicit WARNING when running without real survey data files
"""

import numpy as np
import sys
import os
import warnings

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sct_core import CAR_predictions
from likelihoods.desi_bao import DESI_BAO_Likelihood
from likelihoods.des_y6_lensing import DES_Y6_Likelihood
from likelihoods.hsc_kids_lensing import HSC_KiDS_Likelihood
from likelihoods.plank_cmb import Planck_Likelihood     # file is plank_cmb.py


class CombinedLikelihood:
    """
    Combined likelihood for all CAR datasets.

    When real survey data files are not present in data/, mock data is used.
    Mock data is generated from CAR predictions — results will be circular
    and should NOT be used to claim Bayesian evidence. A warning is printed
    when mock data is active.
    """

    def __init__(self, data_dir=None):
        self.desi     = DESI_BAO_Likelihood(data_dir)
        self.des_y6   = DES_Y6_Likelihood(data_dir)
        self.hsc_kids = HSC_KiDS_Likelihood(data_dir)
        self.planck   = Planck_Likelihood(data_dir)
        self.ndata    = (self.desi.ndata + self.des_y6.ndata
                         + self.hsc_kids.ndata + self.planck.ndata)

        # Warn if any module is using circular mock data
        if any([self.desi.using_mock, self.des_y6.using_mock,
                self.hsc_kids.using_mock, self.planck.using_mock]):
            warnings.warn(
                "\n\n  *** CIRCULAR DATA WARNING ***\n"
                "  One or more likelihood modules are using mock data\n"
                "  generated from CAR predictions. Chi-squared and Bayes\n"
                "  factor results will be circular and cannot be used to\n"
                "  claim evidence in favour of the CAR framework.\n"
                "  Download real survey data files — see data/download_data.sh\n",
                UserWarning, stacklevel=2
            )

    def log_likelihood(self, Omega_m=0.315, R_b=None):
        """
        Combined log-likelihood for CAR framework.

        Parameters
        ----------
        Omega_m : float   Total matter density
        R_b     : float   CAR coherence parameter. Defaults to R_B_DERIVED=0.2545
                          (Paper 17 v4.0 Section 11.6 — derived, not matched).
                          DO NOT pass 0.260 — that was the legacy matched value.

        Returns
        -------
        logL : float   Total log-likelihood
        """
        from sct_core import R_B_DERIVED
        if R_b is None:
            R_b = R_B_DERIVED  # 0.2545 derived (Paper 17 v4.8 Section 11.6)
        preds = CAR_predictions(Omega_m=Omega_m)

        # Extract parameters — use CAMB-verified values
        r_d        = preds['r_d_Mpc']          # 161.4 Mpc (canonical CAR, v4.8.1 audit)
        H0         = preds['H0_km_s_Mpc']       # 70.4 (CAMB)
        S8         = preds['S8']                # 0.783 (analytic, verified)
        IA_bias    = preds['IA_bias']            # 1.0848 (analytic, v4.8.1 audit)
        theta_star = preds['theta_star']         # 1.04105 (100×θ*, not radians)

        logL_desi     = self.desi.log_likelihood(r_d, Omega_m)
        logL_des_y6   = self.des_y6.log_likelihood(S8)
        logL_hsc_kids = self.hsc_kids.log_likelihood(S8)
        logL_planck   = self.planck.log_likelihood(theta_star)

        return logL_desi + logL_des_y6 + logL_hsc_kids + logL_planck

    def compute_chi2(self, Omega_m=0.315, R_b=None):
        return -2.0 * self.log_likelihood(Omega_m, R_b)


def car_loglike(theta):
    """
    CAR log-likelihood for use with PolyChord or other samplers.

    Parameters
    ----------
    theta : array   [Omega_m, R_b]
    """
    Omega_m, R_b = theta
    like = CombinedLikelihood()
    return like.log_likelihood(Omega_m, R_b)


if __name__ == '__main__':
    like   = CombinedLikelihood()
    logL   = like.log_likelihood()
    chi2   = like.compute_chi2()
    print(f'CAR combined log-likelihood : {logL:.2f}')
    print(f'CAR combined chi-squared    : {chi2:.2f}')
    print(f'Total data points           : {like.ndata}')
