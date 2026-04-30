"""
tests/test_sct_core.py
======================
Unit tests for the canonical CAR core calculator (sct_core.py v4.8.1).

These tests verify what the canonical code computes from the canonical
Paper 17 v4.8 §11.6 derived value R_B_DERIVED = 0.2545.

Audited: April 2026 (v4.8.1 NLA recursive audit)
"""

import pytest
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import (
    CAR_predictions, lcdm_reference,
    R_b_of_z, cs_CAR, cs_LCDM, cs_squared, E_of_z,
    compute_S8_analytic, compute_r_d_integral,
    compute_H0_from_theta_and_rd,
    BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2,
    PLANCK_OMEGA_M, PLANCK_Z_STAR, PLANCK_Z_DRAG,
    PLANCK_THETA_STAR, PLANCK_S8,
    R_B_DERIVED, R_B_UNCERTAINTY, R_B_LEGACY_OBS,
    R_D_DERIVED, R_D_UNCERTAINTY,
    N_EFF_SCT, N_EFF_SM, N_EFF_UNCERTAINTY, CMB_S4_SIGMA_NEFF,
    C_HAT_BG_DERIVED,
)


# ── R_b derived constants ──────────────────────────────────────────────────

class TestRbDerivedConstants:
    """Verify Paper 17 v4.8 §11.6 derived constants are correctly set."""

    def test_R_b_derived_value(self):
        assert abs(R_B_DERIVED - 0.2545) < 1e-6

    def test_R_b_uncertainty(self):
        assert abs(R_B_UNCERTAINTY - 0.032) < 1e-6

    def test_R_b_legacy_observed_distinguished(self):
        """The legacy matched value 0.260 is preserved as a reference,
        but the canonical computations all use R_B_DERIVED = 0.2545."""
        assert abs(R_B_LEGACY_OBS - 0.260) < 1e-6
        assert R_B_LEGACY_OBS != R_B_DERIVED

    def test_R_b_agreement_with_observed_within_uncertainty(self):
        """Audit closure: derived 0.2545 ± 0.032 agrees with observed 0.260
        at well below 1σ — closing the v1.0/v2.0 circularity."""
        sigma = abs(R_B_DERIVED - R_B_LEGACY_OBS) / R_B_UNCERTAINTY
        assert sigma < 1.0, f"Derived-vs-observed σ = {sigma:.3f}"


# ── R_b evolution ──────────────────────────────────────────────────────────

class TestRbEvolution:
    def test_R_b_at_z0(self):
        assert abs(R_b_of_z(0.0) - R_B_DERIVED) < 1e-10

    def test_R_b_evolution_form(self):
        """R_b(z) = R_B_DERIVED / (1+z)."""
        for z in [0.0, 0.5, 1.0, 100.0, 1089.0, 1e6]:
            assert abs(R_b_of_z(z) - R_B_DERIVED/(1+z)) < 1e-10

    def test_R_b_at_zstar_essentially_zero(self):
        """At z* ≈ 1089, R_b ≈ 0.000234 — photon limit recovered."""
        Rb_zstar = R_b_of_z(PLANCK_Z_STAR)
        assert Rb_zstar < 1e-3
        assert abs(Rb_zstar - R_B_DERIVED/(1 + PLANCK_Z_STAR)) < 1e-10


# ── Sound speed ────────────────────────────────────────────────────────────

class TestSoundSpeed:
    def test_cs_CAR_at_z0(self):
        """At z=0, c_s²/c² = (1 + 0.2545)/3 = 0.41817."""
        cs2 = cs_CAR(0.0)**2
        assert abs(cs2 - 0.41817) < 1e-4

    def test_cs_CAR_at_z_drag_photon_limit(self):
        """At z_drag ≈ 1060, R_b ≈ 0 → c_s² ≈ 1/3."""
        cs2 = cs_CAR(PLANCK_Z_DRAG)**2
        assert abs(cs2 - 1.0/3.0) < 1e-3

    def test_cs_CAR_high_z_photon_limit(self):
        """At very high z, c_s² → 1/3 (R_b → 0)."""
        cs2 = cs_CAR(1e8)**2
        assert abs(cs2 - 1.0/3.0) < 1e-5

    def test_cs_squared_helper(self):
        """cs_squared(R_b, z) = (1 + R_b/(1+z))/3."""
        for z in [0.0, 1.0, 100.0]:
            for Rb in [0.2545, 0.260, 0.5]:
                expected = (1.0 + Rb/(1+z))/3.0
                assert abs(cs_squared(Rb, z) - expected) < 1e-10


# ── S8 prediction ──────────────────────────────────────────────────────────

class TestS8:
    def test_S8_analytic_formula(self):
        result = compute_S8_analytic()
        expected = PLANCK_S8 * (1 + R_B_DERIVED/3)**(-0.5)
        assert abs(result['S8_analytic'] - expected) < 1e-10

    def test_S8_numerical_offset(self):
        """Numerical S8 = analytic - 0.015 (CAMB Jeans + bary correction)."""
        result = compute_S8_analytic()
        assert abs(result['S8_numeric'] - result['S8_analytic'] - (-0.015)) < 1e-10

    def test_S8_canonical_value(self):
        result = compute_S8_analytic()
        # 0.832 × (1 + 0.2545/3)^(-0.5) = 0.832 × 0.96010 = 0.79881
        assert abs(result['S8_analytic'] - 0.79881) < 1e-4
        assert abs(result['S8_numeric']  - 0.78381) < 1e-4

    def test_S8_uncertainty_propagation(self):
        result = compute_S8_analytic()
        # dS8/dR_b = -0.832/6 × (1 + R_b/3)^(-3/2) × dR_b
        expected_dS8 = (PLANCK_S8 / 6.0) * (1 + R_B_DERIVED/3)**(-1.5) * R_B_UNCERTAINTY
        assert abs(result['dS8_analytic'] - expected_dS8) < 1e-6


# ── b_IA ───────────────────────────────────────────────────────────────────

class TestIABias:
    def test_b_IA_formula(self):
        bIA = 1.0 + R_B_DERIVED / 3.0
        assert abs(bIA - 1.08483) < 1e-4

    def test_b_IA_canonical_constant(self):
        assert abs(C_HAT_BG_DERIVED - 1.08483) < 1e-4

    def test_b_IA_in_predictions(self):
        preds = CAR_predictions()
        assert abs(preds['IA_bias'] - 1.08483) < 1e-4
        assert abs(preds['C_hat_bg'] - 1.08483) < 1e-4


# ── r_d (canonical CAR) ────────────────────────────────────────────────────

class TestSoundHorizon:
    def test_r_d_derived_constant(self):
        """R_D_DERIVED is the canonical CAR r_d (audit v4.8.1 = 161.4 Mpc)."""
        assert abs(R_D_DERIVED - 161.4) < 0.5

    def test_r_d_integral_at_H0_70_4(self):
        """compute_r_d_integral(70.4) should give approximately 161.4 Mpc."""
        rd = compute_r_d_integral(70.4)
        assert abs(rd - 161.4) < 1.0

    def test_r_d_does_not_close_DESI_tension(self):
        """Audit-disclosed honest finding: CAR r_d is NOT 147 Mpc (DESI value)."""
        assert R_D_DERIVED > 155.0, \
            "Canonical CAR r_d should be ≈ 161 Mpc, not the DESI value"


# ── N_eff CMB-S4 prediction ────────────────────────────────────────────────

class TestNeffPrediction:
    def test_N_eff_SCT_value(self):
        assert abs(N_EFF_SCT - 2.514) < 1e-3

    def test_N_eff_SM_value(self):
        assert abs(N_EFF_SM - 3.046) < 1e-3

    def test_N_eff_separation_significant(self):
        """Distance between SCT and SM values, in CMB-S4 σ units."""
        sep = abs(N_EFF_SCT - N_EFF_SM) / CMB_S4_SIGMA_NEFF
        assert sep > 15.0, f"CMB-S4 separation {sep:.1f}σ — should be > 15σ"

    def test_N_eff_signs_opposite(self):
        """SCT predicts N_eff < 3.000 < N_eff_SM — decisive test."""
        assert N_EFF_SCT < 3.0
        assert N_EFF_SM   > 3.0


# ── CAR_predictions integration test ───────────────────────────────────────

class TestCARPredictions:
    @pytest.fixture(scope="class")
    def preds(self):
        return CAR_predictions()

    def test_returns_dict(self, preds):
        assert isinstance(preds, dict)

    def test_all_canonical_keys_present(self, preds):
        for key in ['R_b0', 'cs2', 'C_hat_bg', 'S8', 'S8_analytic',
                    'IA_bias', 'r_d_derived_Mpc', 'N_eff_SCT',
                    'N_eff_CMB_S4_sigma']:
            assert key in preds, f"missing key {key}"

    def test_consistency_with_constants(self, preds):
        """The preds dict should be exactly consistent with the module
        constants — no hidden modifications inside CAR_predictions."""
        assert abs(preds['R_b0']    - R_B_DERIVED) < 1e-10
        assert abs(preds['cs2']     - (1.0 + R_B_DERIVED)/3.0) < 1e-10
        assert abs(preds['C_hat_bg'] - (1.0 + R_B_DERIVED/3.0)) < 1e-10

    def test_S8_analytic_minus_offset_equals_numerical(self, preds):
        assert abs(preds['S8_analytic'] - preds['S8'] - 0.015) < 1e-10


# ── Reproducibility ────────────────────────────────────────────────────────

class TestReproducibility:
    def test_predictions_deterministic(self):
        """Calling CAR_predictions() twice should give identical numbers."""
        p1 = CAR_predictions()
        p2 = CAR_predictions()
        for k in p1:
            v1, v2 = p1[k], p2[k]
            if isinstance(v1, (int, float)) and isinstance(v2, (int, float)):
                assert abs(v1 - v2) < 1e-12, f"{k}: {v1} != {v2}"
            else:
                assert v1 == v2, f"{k}: {v1} != {v2}"

    def test_legacy_R_b_comparison_does_not_corrupt_state(self):
        """Calling with explicit R_b argument should not modify defaults."""
        legacy = CAR_predictions(R_b=R_B_LEGACY_OBS)
        canonical = CAR_predictions()
        assert abs(legacy['R_b0']    - R_B_LEGACY_OBS) < 1e-10
        assert abs(canonical['R_b0'] - R_B_DERIVED)    < 1e-10
        # And calling canonical AGAIN should still be canonical
        canonical2 = CAR_predictions()
        assert abs(canonical2['R_b0'] - R_B_DERIVED) < 1e-10


# ── Audit findings — explicit regression tests ──────────────────────────────

class TestAuditFindings:
    """Regression tests for the v4.8.1 NLA recursive audit findings."""

    def test_finding_1_no_R_b_0_260_in_canonical_path(self):
        """Audit finding #1: R_b = 0.260 must NOT be the canonical value."""
        assert R_B_DERIVED != 0.260

    def test_finding_4_r_d_corrected_from_146_8(self):
        """Audit finding #4: r_d corrected from 146.8 to 161.4 Mpc."""
        assert R_D_DERIVED > 155.0
        assert R_D_DERIVED < 165.0

    def test_finding_7_b_IA_corrected_from_1_087(self):
        """Audit finding #7: b_IA corrected to 1.0848 (was 1.087)."""
        bIA = 1.0 + R_B_DERIVED / 3.0
        assert abs(bIA - 1.0848) < 0.0005


# ── ΛCDM reference (for plotting & comparison) ─────────────────────────────

class TestLCDMReference:
    def test_lcdm_reference_returns_dict(self):
        ref = lcdm_reference()
        for key in ['cs2', 'b_IA', 'S8', 'r_d', 'H0', 'N_eff']:
            assert key in ref, f"missing {key}"

    def test_lcdm_cs2_is_one_third(self):
        assert abs(lcdm_reference()['cs2'] - 1.0/3.0) < 1e-10
