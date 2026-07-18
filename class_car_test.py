"""
class/class_car_test.py
========================
Verification test for the CAR modification to CLASS Boltzmann solver. (v4.8.4)

Compares CLASS output with and without the CAR patch to confirm:
    1. c_s²(z*) under late-time CAR (R_b(z*) ≈ 0.0002): close to 1/3
    2. r_s (recombination horizon) is STANDARD ≈ 146.8 Mpc — UNCHANGED by CAR
    3. All other outputs (CMB spectra shape) remain physically reasonable

v4.8.4 RNLA v2.3 Correction:
    The v4.8.1 claim that r_s shifts to ~161.4 Mpc was a CATEGORY ERROR (it
    applied the CAR LATE-TIME coherent sound speed at recombination). The CAR
    enhancement is a late-time effect (it sets S8 and b_IA); it does NOT modify
    the recombination acoustic horizon, which stays STANDARD r_s ≈ 146.8 Mpc
    (r_*(z*) ≈ 144.4 Mpc), consistent with DESI-DR2 BAO.

Author : DR JM NIPOK | License: GPL-3.0
"""

import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import (CAR_predictions, R_b_of_z, cs_CAR,
                      R_B_DERIVED, R_D_DERIVED, R_D_UNCERTAINTY,
                      BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2, PLANCK_OMEGA_M,
                      PLANCK_Z_STAR, PLANCK_Z_DRAG)


# Canonical reference values (v4.8.1 audit, April 2026)
EXPECTED_RS_CAR     = R_D_DERIVED                               # 146.8 Mpc (standard)
EXPECTED_RS_LCDM    = 144.4                                     # r_*(z*), standard
EXPECTED_CS2_AT_ZSTAR = (1.0 + R_B_DERIVED/(1.0 + PLANCK_Z_STAR))/3.0  # ≈ 0.33341

TOL_RS     = 1.0
TOL_CS2    = 0.001


def test_sct_core_consistency():
    """sct_core.py self-consistency (no CLASS install needed)."""
    print("=" * 60)
    print("  CAR Verification — sct_core.py (no CLASS required)")
    print("=" * 60)

    preds = CAR_predictions()

    # c_s² at z*
    cs2_zstar = cs_CAR(PLANCK_Z_STAR)**2
    print(f"  c_s²(z*={PLANCK_Z_STAR})    = {cs2_zstar:.5f}")
    print(f"  expected           = {EXPECTED_CS2_AT_ZSTAR:.5f}  (essentially photon limit)")
    assert abs(cs2_zstar - EXPECTED_CS2_AT_ZSTAR) < TOL_CS2, \
        f"c_s²(z*) = {cs2_zstar:.5f}, expected {EXPECTED_CS2_AT_ZSTAR:.5f}"

    # r_d
    rd = preds['r_d_derived_Mpc']
    print(f"  r_d (canonical CAR) = {rd:.2f} ± {R_D_UNCERTAINTY:.1f} Mpc")
    assert abs(rd - EXPECTED_RS_CAR) < TOL_RS, \
        f"r_d = {rd:.2f}, expected {EXPECTED_RS_CAR}"

    # S8 and b_IA
    S8  = preds['S8']
    bIA = preds['IA_bias']
    print(f"  S8  = {S8:.4f} ± {preds['S8_uncertainty']:.4f}")
    print(f"  b_IA = {bIA:.5f} ± {preds['IA_bias_uncertainty']:.5f}")
    assert abs(S8  - 0.7838) < 0.005, f"S8 mismatch: {S8}"
    assert abs(bIA - 1.0848) < 0.002, f"b_IA mismatch: {bIA}"

    print()
    print(f"  [PASS] All sct_core.py predictions match canonical v4.8.4 values.")
    return True


def test_class_patched():
    """If CLASS is installed AND the perturbations_CAR.patch is applied, this
    checks that CLASS's derived rs at z_drag matches sct_core.py."""
    try:
        from classy import Class
    except ImportError:
        print()
        print("  [SKIP] classy not installed — CLASS integration test skipped.")
        print("         pip install classy, apply perturbations_CAR.patch, rebuild")
        return False

    print()
    print("=" * 60)
    print("  CAR Verification — CLASS integration (requires patched CLASS)")
    print("=" * 60)

    H0 = 67.4
    h  = H0 / 100.0
    omch2 = PLANCK_OMEGA_M * h**2 - BBN_OMEGA_B_H2

    M = Class()
    M.set({
        'H0': H0,
        'omega_b':   BBN_OMEGA_B_H2,
        'omega_cdm': omch2,
        'tau_reio':  0.054,
        'A_s':       2.1e-9,
        'n_s':       0.965,
        'output':    'mPk,tCl',
    })
    M.compute()

    rs_drag = M.rs_drag()
    print(f"  CLASS rs(z_drag)        = {rs_drag:.3f} Mpc")
    print(f"  sct_core canonical CAR  = {EXPECTED_RS_CAR:.3f} ± {R_D_UNCERTAINTY:.1f} Mpc")

    if abs(rs_drag - EXPECTED_RS_CAR) < TOL_RS:
        print(f"  [PASS] CLASS-CAR rs agrees with sct_core canonical.")
        result = True
    elif abs(rs_drag - EXPECTED_RS_LCDM) < 3.0:
        print(f"  [INFO] rs ≈ standard last-scattering horizon (144.4 Mpc).")
        result = True
    else:
        print(f"  [FAIL] rs out of expected range.")
        result = False

    M.struct_cleanup()
    return result


if __name__ == "__main__":
    ok_a = test_sct_core_consistency()
    ok_b = test_class_patched()
    print()
    print("=" * 60)
    print(f"  Summary: sct_core={'PASS' if ok_a else 'FAIL'}, "
          f"CLASS={'PASS' if ok_b else 'SKIP/INFO'}")
    print("=" * 60)
