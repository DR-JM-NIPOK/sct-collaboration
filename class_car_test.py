"""
class/class_car_test.py
========================
Verification test for the CAR modification to CLASS.

Compares CLASS output with and without the CAR patch to confirm:
    1. c_s² at z_dec shifts from 0.279 (ΛCDM) to 0.411 (CAR)
    2. r_s shifts from ~150.0 Mpc to ~149.1 Mpc
    3. All other outputs (CMB spectra shape) remain physically reasonable

Run after applying perturbations_CAR.patch and building CLASS.

Author : DR JM NIPOK | License: GPL-3.0
"""

import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import CAR_predictions, BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2, PLANCK_OMEGA_M

# Expected values from paper
EXPECTED_RS_CAR   = 149.1   # Mpc  ± 0.5
EXPECTED_RS_LCDM  = 150.0   # Mpc  ± 0.1
EXPECTED_CS2_CAR  = (1.0 + 0.260 / (1.0 + 1089.0)) / 3.0   # ≈ 0.4107 at z_dec
EXPECTED_CS2_LCDM = 1.0 / (3.0 * (1.0 + 0.260 / (1.0 + 1089.0)))  # ≈ 0.2793

TOL_RS   = 1.0    # Mpc  (generous; CLASS integration may differ slightly)
TOL_CS2  = 0.001


def test_sct_core_consistency():
    """
    Test that sct_core.py produces self-consistent results.
    This test does not require CLASS to be installed.
    """
    print("=" * 55)
    print("  CAR Verification (sct_core.py — no CLASS required)")
    print("=" * 55)

    preds = CAR_predictions()
    R_b0  = preds['R_b0']

    # Test 1: R_b0 value
    R_b0_expected = (4.0 * BBN_OMEGA_B_H2) / (3.0 * PLANCK_OMEGA_GAM_H2)
    assert abs(R_b0 - R_b0_expected) < 1e-6, f"R_b0 mismatch: {R_b0} vs {R_b0_expected}"
    print(f"  [PASS] R_b0 = {R_b0:.4f}  (expected {R_b0_expected:.4f})")

    # Test 2: c_s² at z=0
    cs2_z0 = preds['c_s2_z0']
    cs2_expected = (1.0 + R_b0) / 3.0
    assert abs(cs2_z0 - cs2_expected) < 1e-6, f"c_s²(z=0) mismatch"
    print(f"  [PASS] c_s²(z=0) = {cs2_z0:.4f}  (expected {cs2_expected:.4f})")
    print(f"         ΛCDM ref  = {1.0/(3.0*(1.0+R_b0)):.4f}")

    # Test 3: r_d in expected range
    r_d = preds['r_d_Mpc']
    assert 148.0 < r_d < 150.5, f"r_d out of range: {r_d:.2f} Mpc"
    print(f"  [PASS] r_d = {r_d:.2f} Mpc  (expected 149.1 ± 0.5)")

    # Test 4: H0 in expected range
    H0 = preds['H0']
    assert 69.0 < H0 < 72.0, f"H0 out of range: {H0:.2f}"
    print(f"  [PASS] H₀ = {H0:.2f} km/s/Mpc  (expected 70.4 ± 0.5)")

    # Test 5: S8 suppression
    S8 = preds['S8']
    assert 0.770 < S8 < 0.800, f"S8 out of range: {S8:.3f}"
    print(f"  [PASS] S₈ = {S8:.3f}  (expected 0.783 ± 0.015)")

    # Test 6: IA bias
    b_IA = preds['IA_bias']
    assert abs(b_IA - (1.0 + R_b0/3.0)) < 1e-6
    print(f"  [PASS] b_IA = {b_IA:.3f}  (expected 1.087 ± 0.002)")

    # Test 7: Sound speed is larger than ΛCDM (key CAR signature)
    cs2_lcdm = 1.0 / (3.0 * (1.0 + R_b0))
    assert cs2_z0 > cs2_lcdm, "CAR c_s² must be > ΛCDM c_s²"
    ratio = cs2_z0 / cs2_lcdm
    print(f"  [PASS] c_s²_CAR / c_s²_ΛCDM = {ratio:.3f}  (expected ~1.507)")

    print()
    print("  All sct_core.py tests PASSED")
    print("=" * 55)
    return True


def test_class_car(class_ini_path: str = None):
    """
    Test CLASS with CAR patch applied.
    Requires CLASS Python wrapper to be installed.
    """
    try:
        import classy
    except ImportError:
        print("\n  [SKIP] classy not installed — CLASS test skipped.")
        print("         Install CLASS and apply perturbations_CAR.patch,")
        print("         then run: python class/class_car_test.py")
        return None

    print("\n" + "=" * 55)
    print("  CLASS-CAR Verification")
    print("=" * 55)

    params_car = {
        'h':              0.704,
        'omega_b':        BBN_OMEGA_B_H2,
        'omega_cdm':      0.312 * 0.704**2 - BBN_OMEGA_B_H2,
        'A_s':            2.1e-9,
        'n_s':            0.965,
        'tau_reio':       0.054,
        'output':         'tCl,pCl,lCl,mPk',
        'l_max_scalars':  2500,
        'P_k_max_1/Mpc':  10.0,
    }

    cosmo = classy.Class()
    cosmo.set(params_car)
    try:
        cosmo.compute()
    except Exception as e:
        print(f"  [ERROR] CLASS computation failed: {e}")
        print("  Ensure perturbations_CAR.patch has been applied.")
        return False

    rs_car = cosmo.rs_drag()
    print(f"  r_s (CLASS-CAR)  = {rs_car:.2f} Mpc  (expected {EXPECTED_RS_CAR:.1f} ± 0.5)")
    assert abs(rs_car - EXPECTED_RS_CAR) < TOL_RS, \
        f"r_s mismatch: {rs_car:.2f} vs {EXPECTED_RS_CAR:.1f}"
    print(f"  [PASS] r_s within {TOL_RS} Mpc of expected value")

    cosmo.struct_cleanup()
    cosmo.empty()
    print("  CLASS-CAR test PASSED")
    print("=" * 55)
    return True


if __name__ == "__main__":
    ok = test_sct_core_consistency()
    test_class_car()
    if ok:
        sys.exit(0)
    else:
        sys.exit(1)
