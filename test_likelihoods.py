"""
tests/test_likelihoods.py
==========================
Unit tests for all CAR likelihood modules.
No external data files required — tests use mock data.

Author : DR JM NIPOK | License: GPL-3.0
"""

import pytest
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import CAR_predictions


@pytest.fixture(scope="module")
def params():
    return CAR_predictions()


class TestDESIDR2:
    def test_import(self):
        from likelihoods.desi_dr2_bao import DESIDR2BAOLikelihood
        lik = DESIDR2BAOLikelihood(verbose=False)
        assert lik.name == "DESI-DR2 BAO"

    def test_log_like_finite(self, params):
        from likelihoods.desi_dr2_bao import DESIDR2BAOLikelihood
        lik = DESIDR2BAOLikelihood(verbose=False)
        lnL = lik.log_like(params)
        assert np.isfinite(lnL), f"log-like = {lnL}"

    def test_chi2_per_dof_reasonable(self, params):
        from likelihoods.desi_dr2_bao import DESIDR2BAOLikelihood
        lik = DESIDR2BAOLikelihood(verbose=False)
        c = lik.chi2_per_dof(params)
        # v4.8.1 audit: canonical CAR r_d ≈ 161.4 vs DESI ≈ 147 — chi2/dof
        # is large (CAR doesn't fit DESI). The honest test is that the
        # number is finite and computable, not that it's small.
        assert np.isfinite(c) and c > 0, f"chi2/dof = {c} should be positive finite"

    def test_tension_summary(self, params):
        from likelihoods.desi_dr2_bao import DESIDR2BAOLikelihood
        lik = DESIDR2BAOLikelihood()
        t = lik.tension_summary(params)
        assert 'tension_sigma' in t
        # v4.8.1 audit: CAR r_d ≈ 161.4 (NOT 149.1). DESI = 147.0 ± 1.0
        # → tension ≈ 14σ. The audit-disclosed honest finding is that
        # CAR does NOT close the DESI-DR2 BAO tension. This assertion
        # documents that fact for future regressions.
        assert t['tension_sigma'] > 10.0, \
            f"CAR r_d should be in HIGH tension with DESI; got {t['tension_sigma']:.1f}σ"
        assert t['tension_sigma'] < 20.0, \
            f"Tension > 20σ would indicate computation drift; got {t['tension_sigma']:.1f}σ"


class TestDESY6:
    def test_import(self):
        from likelihoods.des_y6_3x2pt import DESY6ThreeTwoPointLikelihood
        lik = DESY6ThreeTwoPointLikelihood(verbose=False)
        assert lik.n_data == 1800

    def test_log_like_finite(self, params):
        from likelihoods.des_y6_3x2pt import DESY6ThreeTwoPointLikelihood
        lik = DESY6ThreeTwoPointLikelihood(verbose=False)
        lnL = lik.log_like(params)
        assert np.isfinite(lnL)

    def test_S8_summary(self, params):
        from likelihoods.des_y6_3x2pt import DESY6ThreeTwoPointLikelihood
        lik = DESY6ThreeTwoPointLikelihood()
        s = lik.S8_posterior_summary(params)
        assert s['tension_sigma'] < 0.5   # CAR matches DES-Y6 to <0.5σ


class TestHSCY3:
    def test_import(self):
        from likelihoods.hsc_y3_wl import HSCY3WeakLensingLikelihood
        lik = HSCY3WeakLensingLikelihood()
        assert lik.n_data == 576

    def test_log_like_finite(self, params):
        from likelihoods.hsc_y3_wl import HSCY3WeakLensingLikelihood
        lik = HSCY3WeakLensingLikelihood(verbose=False)
        lnL = lik.log_like(params)
        assert np.isfinite(lnL)


class TestKiDSDR5:
    def test_import(self):
        from likelihoods.kids_dr5_wl import KiDSDR5WeakLensingLikelihood
        lik = KiDSDR5WeakLensingLikelihood()
        assert lik.n_data == 912

    def test_log_like_finite(self, params):
        from likelihoods.kids_dr5_wl import KiDSDR5WeakLensingLikelihood
        lik = KiDSDR5WeakLensingLikelihood(verbose=False)
        lnL = lik.log_like(params)
        assert np.isfinite(lnL)

    def test_S8_ratio(self):
        from likelihoods.kids_dr5_wl import KiDSDR5WeakLensingLikelihood
        lik = KiDSDR5WeakLensingLikelihood()
        r = lik.S8_ratio_test()
        # Observed ratio ~0.947, CAR prediction ~0.959, within ~0.7σ
        assert r['tension_sigma'] < 2.0


class TestPlanckPR4:
    def test_import(self):
        from likelihoods.planck_pr4_cmb import PlanckPR4CMBLikelihood
        lik = PlanckPR4CMBLikelihood()
        assert "Planck" in lik.name

    def test_log_like_finite(self, params):
        from likelihoods.planck_pr4_cmb import PlanckPR4CMBLikelihood
        lik = PlanckPR4CMBLikelihood(verbose=False)
        lnL = lik.log_like(params)
        assert np.isfinite(lnL)


class TestCombinedLikelihood:
    def test_combined_log_like(self, params):
        """Combined log-likelihood should be finite and negative."""
        from likelihoods import (DESIDR2BAOLikelihood, DESY6ThreeTwoPointLikelihood,
                                  HSCY3WeakLensingLikelihood, KiDSDR5WeakLensingLikelihood,
                                  PlanckPR4CMBLikelihood)
        liks = [
            DESIDR2BAOLikelihood(),
            DESY6ThreeTwoPointLikelihood(),
            HSCY3WeakLensingLikelihood(),
            KiDSDR5WeakLensingLikelihood(),
            PlanckPR4CMBLikelihood(),
        ]
        total_lnL = sum(lik.log_like(params) for lik in liks)
        assert np.isfinite(total_lnL)
        assert total_lnL < 0


"""
tests/test_sound_horizon.py
============================
Tests for sound horizon numerical integration.
"""


class TestSoundHorizonIntegration:
    def test_integrand_positive(self):
        from sct_core import cs_squared, hubble_factor, R_B_DERIVED
        # v4.8.1 audit: use canonical derived R_b, not the broken
        # 4·Ω_b/(3·Ω_γ) formula (= 1196.9, never the SCT R_b).
        R_b0 = R_B_DERIVED
        for z in [100, 500, 1000, 1089, 2000]:
            cs2 = cs_squared(R_b0, z)
            Hz  = hubble_factor(z, 0.312, 9e-5)
            integrand = np.sqrt(cs2) / Hz
            assert integrand > 0, f"Integrand non-positive at z={z}"

    def test_integral_convergence(self):
        from scipy.integrate import quad
        from sct_core import cs_squared, hubble_factor, R_B_DERIVED
        # v4.8.1 audit: canonical R_b
        R_b0 = R_B_DERIVED
        def integrand(z):
            return np.sqrt(cs_squared(R_b0, z)) / hubble_factor(z, 0.312, 9e-5)
        # The integrand falls as 1/z² at high z; need a sufficiently high
        # upper limit (~1e7) AND a high comparison limit to see convergence.
        result_1, _ = quad(integrand, 1060, 1e7, limit=400, epsrel=1e-10)
        result_2, _ = quad(integrand, 1060, 1e9, limit=400, epsrel=1e-10)
        rel = abs(result_2 - result_1) / result_1
        assert rel < 0.005, \
            f"Convergence: {rel:.4%} > 0.5%"

    def test_rd_monotone_in_Rb(self):
        from sct_core import compute_sound_horizon, PLANCK_OMEGA_M, PLANCK_Z_STAR
        R_b_values = [0.24, 0.26, 0.28]
        r_d_values = [compute_sound_horizon(R_b, PLANCK_OMEGA_M, 9e-5, PLANCK_Z_STAR)
                      for R_b in R_b_values]
        # Higher R_b → larger c_s → larger r_d
        assert r_d_values[0] < r_d_values[1] < r_d_values[2]
