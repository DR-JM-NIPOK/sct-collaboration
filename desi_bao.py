"""
desi_bao.py — DESI-DR2 BAO likelihood for CAR framework
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - Mock data fallback now clearly labelled as CIRCULAR
  - using_mock flag added so CombinedLikelihood can warn users
  - Real data path documented with download instructions
  - Mock data uses ΛCDM fiducial (r_d=147.1) not CAR prediction
    to avoid circular chi-squared = 0 results
"""

import numpy as np
import os
from scipy.integrate import quad


class DESI_BAO_Likelihood:
    """
    DESI-DR2 BAO likelihood.

    Loads real data from data/desi_dr2_bao.txt if present.
    Falls back to a ΛCDM-fiducial mock ONLY for code testing.
    Do NOT use mock data to compute Bayesian evidence.

    Real data download:
        bash data/download_data.sh
    """

    def __init__(self, data_dir=None):
        if data_dir is None:
            self.data_dir = os.path.join(
                os.path.dirname(__file__), '..', 'data')
        else:
            self.data_dir = data_dir

        self.ndata      = 32
        self.using_mock = False
        self.load_data()

    def load_data(self):
        data_file = os.path.join(self.data_dir, 'desi_dr2_bao.txt')
        cov_file  = os.path.join(self.data_dir, 'desi_dr2_cov.txt')

        if os.path.exists(data_file) and os.path.exists(cov_file):
            data     = np.loadtxt(data_file)
            self.z             = data[:, 0]
            self.D_V_over_r_d  = data[:, 1]
            self.sigma         = data[:, 2] if data.shape[1] > 2 else None
            self.cov           = np.loadtxt(cov_file)
            self.using_mock    = False
        else:
            # ── MOCK DATA — FOR CODE TESTING ONLY ─────────────────────────
            # Uses ΛCDM fiducial r_d = 147.1 Mpc (NOT CAR's 149.2 Mpc).
            # Using CAR's own r_d here would make chi-squared ≈ 0 by
            # construction and produce circular Bayesian evidence.
            import warnings
            warnings.warn(
                "DESI: using ΛCDM-fiducial mock data (r_d=147.1 Mpc). "
                "Run data/download_data.sh to get real DESI-DR2 data.",
                UserWarning, stacklevel=2
            )
            self.using_mock = True

            self.z = np.array([
                0.295, 0.510, 0.706, 0.930, 1.317, 1.491, 2.330
            ])

            # ΛCDM fiducial D_V/r_d values (Omega_m=0.315, r_d=147.1 Mpc)
            Omega_m_fid = 0.315
            r_d_fid     = 147.1   # ΛCDM, NOT CAR

            def H_over_H0(z):
                return np.sqrt(Omega_m_fid * (1+z)**3 + (1-Omega_m_fid))

            def D_V_rd(z):
                I, _ = quad(lambda zp: 1.0/H_over_H0(zp), 0, z)
                d_A  = I / (1+z)
                D_V  = (z * d_A**2 / H_over_H0(z))**(1/3)
                return D_V * 2997.9 / r_d_fid   # c/100 Mpc/h → Mpc

            self.D_V_over_r_d = np.array([D_V_rd(zi) for zi in self.z])

            # Published DESI-DR2 uncertainties (approximate)
            sigma_frac = np.array([0.011, 0.008, 0.007, 0.006, 0.006, 0.008, 0.012])
            self.cov   = np.diag((sigma_frac * self.D_V_over_r_d)**2)

    def _model(self, z, r_d, Omega_m):
        """Compute D_V/r_d prediction for given parameters."""
        def H_over_H0(zp):
            return np.sqrt(Omega_m * (1+zp)**3 + (1-Omega_m))

        I, _ = quad(lambda zp: 1.0/H_over_H0(zp), 0, z)
        d_A  = I / (1+z)
        D_V  = (z * d_A**2 / H_over_H0(z))**(1/3)
        return D_V * 2997.9 / r_d

    def compute_chi2(self, r_d, Omega_m):
        model     = np.array([self._model(zi, r_d, Omega_m) for zi in self.z])
        residuals = self.D_V_over_r_d - model
        inv_cov   = np.linalg.inv(self.cov)
        return float(residuals @ inv_cov @ residuals)

    def log_likelihood(self, r_d, Omega_m):
        return -0.5 * self.compute_chi2(r_d, Omega_m)
