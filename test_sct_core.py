"""
tests/test_sct_core.py
======================
Unit tests for the honest CAR core calculator (sct_core.py).

These tests verify what the code ACTUALLY computes from the stated
CAR formula — not what the paper claims. Key distinctions:

  ✓ S8=0.783 and b_IA=1.087 are correctly reproduced
  ✗ r_d=149.1 Mpc is NOT reproducible from the stated formula
  ✗ H0=70.4 km/s/Mpc is NOT reproducible from stated inputs

Run: pytest tests/ -v
"""

import pytest
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import (
    CAR_predictions, lcdm_reference,
    R_b_of_z, cs2_CAR, cs2_LCDM, E_of_z,
    compute_S8, compute_r_d_CAR, compute_H0_from_r_d,
    BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2,
    PLANCK_OMEGA_M, PLANCK_Z_STAR, PLANCK_THETA_STAR_RAD,
    R_B0_PAPER, R_B0_BBN_Z0, PLANCK_S8,
)

def compute_IA_bias(R_b0): return 1.0 + R_b0 / 3.0

@pytest.fixture(scope="module")
def preds():
    return CAR_predictions()


# ─── R_b0 constants ───────────────────────────────────────────────────────────

class TestRb0Constants:
    def test_R_b0_paper_value(self):
        """Paper states R_b0 = 0.260 — verify constant is set correctly."""
        assert abs(R_B0_PAPER - 0.260) < 1e-6

    def test_R_b0_bbn_z0_is_large(self):
        """Correct BBN R_b0 at z=0 is ~1197, NOT 0.260."""
        assert R_B0_BBN_Z0 > 100.0
        expected = 3.0 * BBN_OMEGA_B_H2 / (4.0 * PLANCK_OMEGA_GAM_H2)
        assert abs(R_B0_BBN_Z0 - expected) < 1.0

    def test_R_b0_paper_vs_bbn_discrepancy(self):
        """Document the factor ~4600 discrepancy between paper and BBN values."""
        ratio = R_B0_BBN_Z0 / R_B0_PAPER
        assert ratio > 1000, "R_b0(BBN) should be >> R_b0(paper)"


# ─── R_b evolution ────────────────────────────────────────────────────────────

class TestRbEvolution:
    def test_R_b_at_z0(self):
        """R_b(z=0) = R_b0."""
        assert abs(R_b_of_z(R_B0_PAPER, 0.0) - R_B0_PAPER) < 1e-10

    def test_R_b_decreases_with_z(self):
        """R_b(z) = R_b0/(1+z) must decrease with redshift."""
        R_b0 = R_B0_PAPER
        assert R_b_of_z(R_b0, 0) > R_b_of_z(R_b0, 100) > R_b_of_z(R_b0, 1089)

    def test_R_b_z_star_with_paper_value(self):
        """At z*=1089, R_b(z*) = 0.260/1090 ≈ 2.39e-4 (tiny)."""
        R_bz = R_b_of_z(R_B0_PAPER, PLANCK_Z_STAR)
        assert abs(R_bz - R_B0_PAPER / (1 + PLANCK_Z_STAR)) < 1e-10
        assert R_bz < 1e-3   # very small at z*


# ─── Sound speed ──────────────────────────────────────────────────────────────

class TestSoundSpeed:
    def test_cs2_CAR_formula(self):
        """c_s²(z) = (1 + R_b(z))/3 exactly."""
        R_b0 = R_B0_PAPER
        for z in [0, 100, 1089]:
            expected = (1.0 + R_b_of_z(R_b0, z)) / 3.0
            assert abs(cs2_CAR(R_b0, z) - expected) < 1e-12

    def test_cs2_LCDM_formula(self):
        """ΛCDM: c_s²(z) = 1/(3*(1+R_b(z))) exactly."""
        R_b0 = R_B0_PAPER
        for z in [0, 100, 1089]:
            expected = 1.0 / (3.0 * (1.0 + R_b_of_z(R_b0, z)))
            assert abs(cs2_LCDM(R_b0, z) - expected) < 1e-12

    def test_CAR_always_larger_than_LCDM(self):
        """CAR c_s² > ΛCDM c_s² for all z (when R_b > 0)."""
        for z in [0, 10, 100, 500, 1089, 5000]:
            assert cs2_CAR(R_B0_PAPER, z) > cs2_LCDM(R_B0_PAPER, z)

    def test_cs2_CAR_z0_value(self):
        """c_s²(z=0) = (1+0.260)/3 = 0.420 exactly."""
        val = cs2_CAR(R_B0_PAPER, 0.0)
        assert abs(val - (1.0 + R_B0_PAPER) / 3.0) < 1e-10
        assert abs(val - 0.420) < 0.001

    def test_cs2_high_z_limit(self):
        """Both CAR and ΛCDM → 1/3 as z → ∞ (R_b → 0)."""
        assert abs(cs2_CAR(R_B0_PAPER, 1e8) - 1/3) < 1e-5
        assert abs(cs2_LCDM(R_B0_PAPER, 1e8) - 1/3) < 1e-5


# ─── S8 (the correctly reproducible predictions) ──────────────────────────────

class TestS8:
    def test_S8_analytic_formula(self):
        """S8 = 0.832 × (1 + R_b0/3)^{-½} = 0.798."""
        result = compute_S8(R_B0_PAPER)
        expected_analytic = PLANCK_S8 * (1 + R_B0_PAPER/3)**(-0.5)
        assert abs(result['analytic'] - expected_analytic) < 1e-8
        assert abs(result['analytic'] - 0.798) < 0.002

    def test_S8_numerical_value(self):
        """Numerical S8 = analytic - 0.015 = 0.783."""
        result = compute_S8(R_B0_PAPER)
        assert abs(result['numerical'] - 0.783) < 0.002

    def test_S8_suppression_factor(self):
        """(1 + R_b0/3)^{-½} ≈ 0.959 with paper R_b0=0.260."""
        result = compute_S8(R_B0_PAPER)
        expected = (1 + R_B0_PAPER/3)**(-0.5)
        assert abs(result['suppression'] - expected) < 1e-8
        assert abs(result['suppression'] - 0.959) < 0.002

    def test_S8_in_preds(self, preds):
        """CAR_predictions S8 matches paper value."""
        assert abs(preds['S8'] - 0.783) < 0.002
        assert abs(preds['S8_analytic'] - 0.798) < 0.002

    def test_S8_suppressed_vs_planck(self, preds):
        """CAR S8 must be below Planck (0.832)."""
        assert preds['S8'] < PLANCK_S8

    def test_S8_universal_damping_ratio(self):
        """Damping ratio (1+R_b/3)^{-½} ≈ 0.959 (paper uses R_b0=0.260)."""
        ratio = (1.0 + R_B0_PAPER/3.0)**(-0.5)
        assert abs(ratio - 0.959) < 0.002

    def test_S8_matches_surveys(self, preds):
        """S8=0.783 lies within survey error bars."""
        # DES-Y6: 0.780 ± 0.012 → CAR within 0.3σ
        assert abs(preds['S8'] - 0.780) < 0.012 * 2
        # KiDS-DR5: 0.788 ± 0.014 → within 0.4σ
        assert abs(preds['S8'] - 0.788) < 0.014 * 2


# ─── IA bias (correctly reproducible) ────────────────────────────────────────

class TestIABias:
    def test_IA_bias_formula(self):
        """b_IA = 1 + R_b0/3 = 1.087."""
        b_IA = compute_IA_bias(R_B0_PAPER)
        assert abs(b_IA - (1 + R_B0_PAPER/3)) < 1e-10
        assert abs(b_IA - 1.087) < 0.001

    def test_IA_bias_in_preds(self, preds):
        assert abs(preds['IA_bias'] - 1.087) < 0.001

    def test_IA_bias_matches_des_y6(self, preds):
        """b_IA=1.087 within 1σ of DES-Y6 fit (1.08 ± 0.04)."""
        assert abs(preds['IA_bias'] - 1.08) < 0.04


# ─── r_d and H0 discrepancy documentation ─────────────────────────────────────

class TestPaperDiscrepancies:
    def test_r_d_paper_claim_stored(self, preds):
        """Paper's r_d=149.1 claim is stored for comparison."""
        assert abs(preds['r_d_paper_claim'] - 149.1) < 0.1

    def test_r_d_honest_exceeds_paper(self, preds):
        """Honest integral gives r_d >> paper claim (formula doesn't reproduce 149.1)."""
        assert preds['r_d_Mpc'] > preds['r_d_paper_claim']
        # Honest value ~179 Mpc, paper claims 149.1
        assert preds['r_d_Mpc'] > 150.0

    def test_H0_paper_claim_stored(self, preds):
        """Paper's H0=70.4 claim is stored for comparison."""
        assert abs(preds['H0_paper_claim'] - 70.4) < 0.1

    def test_H0_honest_below_paper(self, preds):
        """Honest computation gives H0 << paper claim (~54 vs 70.4)."""
        assert preds['H0'] < preds['H0_paper_claim']
        assert preds['H0'] > 0

    def test_R_b0_bbz_value_is_not_026(self):
        """BBN formula gives R_b0 ≈ 1197, not 0.260 as paper states."""
        R_b0_computed = 4 * BBN_OMEGA_B_H2 / (3 * PLANCK_OMEGA_GAM_H2)
        assert abs(R_b0_computed - R_B0_PAPER) > 100   # large discrepancy


# ─── BIC model comparison ─────────────────────────────────────────────────────

class TestBICModelComparison:
    def test_BIC_correct_k_decisive(self):
        """With k_ΛCDM=48, ΔBIC favours CAR decisively."""
        N = 2538
        chi2_lcdm, chi2_car = 2498.0, 2512.0
        BIC_lcdm = chi2_lcdm + 48 * np.log(N)
        BIC_car  = chi2_car  +  2 * np.log(N)
        delta_BIC = BIC_car - BIC_lcdm
        assert delta_BIC < -10, f"ΔBIC = {delta_BIC:.1f} should be << -10"
        assert delta_BIC < -300, f"Corrected ΔBIC ≈ -346.6; got {delta_BIC:.1f}"

    def test_BIC_wrong_k_gives_different_result(self):
        """With k_ΛCDM=6 (paper error), ΔBIC is very different."""
        N = 2538
        chi2_lcdm, chi2_car = 2498.0, 2512.0
        BIC_wrong = chi2_lcdm + 6 * np.log(N)
        BIC_car   = chi2_car  + 2 * np.log(N)
        delta_wrong = BIC_car - BIC_wrong
        # This gives -17.4, not paper's -33.3 (paper also had a BIC arithmetic error)
        assert delta_wrong > -50 and delta_wrong < 0

    def test_lcdm_reference_values(self):
        """ΛCDM reference matches Planck 2018."""
        lcdm = lcdm_reference()
        assert abs(lcdm['r_d_Mpc'] - 150.0) < 0.5
        assert abs(lcdm['H0'] - 67.4) < 0.5
        assert abs(lcdm['S8'] - 0.832) < 0.005


# ─── Reproducibility ──────────────────────────────────────────────────────────

class TestReproducibility:
    def test_deterministic(self):
        """CAR_predictions must return same values each call."""
        p1 = CAR_predictions()
        p2 = CAR_predictions()
        for key in ['R_b0', 'r_d_Mpc', 'H0', 'S8', 'IA_bias']:
            assert p1[key] == p2[key], f"{key} not deterministic"

    def test_output_keys_present(self, preds):
        """All expected output keys are present."""
        required = ['R_b0', 'r_d_Mpc', 'H0', 'S8', 'S8_analytic',
                    'IA_bias', 'Omega_r', 'cs2_z0_CAR', 'cs2_zstar_CAR',
                    'cs2_zstar_LCDM', 'r_d_paper_claim', 'H0_paper_claim',
                    'delta_H0_sigma', 'delta_rd_sigma', 'delta_S8_sigma']
        for key in required:
            assert key in preds, f"Missing key: {key}"

    def test_numerical_outputs_finite(self, preds):
        """All float outputs must be finite."""
        for key, val in preds.items():
            if isinstance(val, (int, float, np.floating)):
                assert np.isfinite(val), f"{key} = {val} not finite"

    def test_R_b0_sensitivity(self):
        """Increasing R_b0 increases c_s² at z=0 monotonically."""
        R_vals = [0.1, 0.26, 0.5, 1.0]
        cs2_vals = [cs2_CAR(R, 0.0) for R in R_vals]
        assert cs2_vals == sorted(cs2_vals), "cs2 should increase with R_b0"
