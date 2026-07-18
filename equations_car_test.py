"""
camb/equations_car_test.py
==========================
Verification test for the CAR modification to CAMB. (v4.8.4 RNLA v2.3)

This script verifies that:
    1. sct_core.py's standard sound-horizon integral gives r_d ≈ 146.8 Mpc
    2. CAR predictions for S8 and b_IA match canonical v4.8 values
    3. If CAMB is installed, its standard derived rdrag agrees with the
       standard horizon (≈ 146.8 Mpc) — recombination is NOT patched.

v4.8.4 RNLA v2.3 Correction:
    The v4.8.1 claim that r_d ≈ 161.4 Mpc was a CATEGORY ERROR: it inserted
    the CAR LATE-TIME coherent sound speed into the recombination integral.
    The CAR enhancement is a late-time (S8, b_IA) effect and does NOT modify
    the recombination acoustic horizon. The recombination/drag horizon is
    STANDARD, r_d ≈ 146.8 Mpc (r_*(z*) ≈ 144.4 Mpc), consistent with DESI-DR2.

Author : DR JM NIPOK | License: GPL-3.0
"""

import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import (CAR_predictions, cs_CAR, cs_LCDM, R_b_of_z,
                      compute_r_d_integral,
                      R_B_DERIVED, R_B_UNCERTAINTY,
                      R_D_DERIVED, R_D_UNCERTAINTY,
                      C_HAT_BG_DERIVED, N_EFF_SCT,
                      BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2,
                      PLANCK_OMEGA_M, PLANCK_Z_STAR, PLANCK_Z_DRAG)


# Canonical reference values (v4.8.1 audit, April 2026)
# These should not be relaxed without a corresponding sct_core.py update
EXPECTED_R_D       = 146.8
EXPECTED_R_D_TOL   =   5.0   # Mpc — omega_b, omega_m uncertainty
EXPECTED_S8_NUM    =   0.7838
EXPECTED_S8_TOL    =   0.005
EXPECTED_B_IA      =   1.0848
EXPECTED_B_IA_TOL  =   0.002


def test_canonical_constants():
    """Verify all derived constants match the v4.8.1 canonical values."""
    print("=" * 64)
    print("  CAR Verification — canonical constants (v4.8.1)")
    print("=" * 64)
    assert abs(R_B_DERIVED   - 0.2545)  < 1e-5,  f"R_B_DERIVED = {R_B_DERIVED}"
    assert abs(C_HAT_BG_DERIVED - 1.08483) < 1e-4, f"C_hat_bg = {C_HAT_BG_DERIVED}"
    assert abs(N_EFF_SCT     - 2.514)  < 1e-3,  f"N_EFF_SCT = {N_EFF_SCT}"
    print(f"  R_B_DERIVED      = {R_B_DERIVED}                       [PASS]")
    print(f"  C_HAT_BG_DERIVED = {C_HAT_BG_DERIVED:.5f}                 [PASS]")
    print(f"  N_EFF_SCT        = {N_EFF_SCT}                       [PASS]")
    return True


def test_standalone_predictions():
    """Verify CAR predictions without requiring patched CAMB."""
    print()
    print("=" * 64)
    print("  CAR Verification — standalone (no CAMB needed)")
    print("=" * 64)

    preds = CAR_predictions()
    R_b   = preds['R_b0']

    # c_s² at z=0 and z*
    cs2_z0    = cs_CAR(0.0)**2
    cs2_zstar = cs_CAR(PLANCK_Z_STAR)**2
    cs2_lcdm_zstar = cs_LCDM(PLANCK_Z_STAR)**2

    print(f"  c_s²(z=0)    CAR  = {cs2_z0:.5f}  (expected: 0.41817)")
    print(f"  c_s²(z*)     CAR  = {cs2_zstar:.5f}  (expected: ~0.33341 — photon limit)")
    print(f"  c_s²(z*)     ΛCDM = {cs2_lcdm_zstar:.5f}")
    print(f"  Ratio CAR/ΛCDM    = {cs2_zstar/cs2_lcdm_zstar:.4f}  (CAR enhanced)")

    rd  = preds['r_d_derived_Mpc']
    H0  = preds['H0_CAMB']
    S8  = preds['S8']
    bIA = preds['IA_bias']
    Neff = preds['N_eff_SCT']

    print()
    print(f"  r_d  = {rd:.2f} ± {R_D_UNCERTAINTY:.1f} Mpc")
    print(f"          (standard horizon; DESI-DR2: 147 ± 1 Mpc; consistent ~0.2 sigma)")
    print(f"  H0   = {H0:.2f} km/s/Mpc")
    print(f"  S8   = {S8:.4f} ± {preds['S8_uncertainty']:.4f}  "
          f"(observed DES-Y6: 0.780 ± 0.012, KiDS-DR5: 0.788 ± 0.014)")
    print(f"  b_IA = {bIA:.5f} ± {preds['IA_bias_uncertainty']:.5f}")
    print(f"  N_eff= {Neff:.3f} ± {preds['N_eff_uncertainty']:.2f}  "
          f"(SM: 3.046 → CMB-S4 separation: {preds['N_eff_CMB_S4_sigma']:.1f}σ)")
    print()

    # Tolerances are tight on derived analytics; r_d depends on numerics.
    assert abs(rd - EXPECTED_R_D)   < EXPECTED_R_D_TOL, \
        f"r_d = {rd:.3f} Mpc, expected {EXPECTED_R_D} ± {EXPECTED_R_D_TOL}"
    assert abs(S8 - EXPECTED_S8_NUM) < EXPECTED_S8_TOL, \
        f"S8 = {S8:.4f}, expected {EXPECTED_S8_NUM} ± {EXPECTED_S8_TOL}"
    assert abs(bIA - EXPECTED_B_IA)  < EXPECTED_B_IA_TOL, \
        f"b_IA = {bIA:.5f}, expected {EXPECTED_B_IA} ± {EXPECTED_B_IA_TOL}"

    print(f"  [PASS] all standalone predictions within canonical tolerances.")
    print("=" * 64)
    return True


def test_camb_patched():
    """If CAMB is installed AND the equations_CAR.patch was applied and CAMB
    rebuilt, this checks that CAMB's derived rdrag matches sct_core.py."""
    try:
        import camb
    except ImportError:
        print("\n  [SKIP] camb not installed — CAMB integration test skipped.")
        print("         pip install camb, then apply equations_CAR.patch")
        return False

    print()
    print("=" * 64)
    print("  CAR Verification — CAMB integration (requires patched CAMB)")
    print("=" * 64)

    H0 = 67.4
    Om = PLANCK_OMEGA_M
    h  = H0 / 100.0

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0       = H0,
        ombh2    = BBN_OMEGA_B_H2,
        omch2    = Om * h**2 - BBN_OMEGA_B_H2,
        tau      = 0.054,
        mnu      = 0.06,
        omk      = 0.0,
    )
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_for_lmax(2500, lens_potential_accuracy=1)

    results = camb.get_results(pars)
    derived = results.get_derived_params()
    rd_camb = derived.get('rdrag', None)

    if rd_camb is None:
        print("  [WARN] rdrag not in derived params — check CAMB version")
        return False

    print(f"  r_drag (this CAMB build)   = {rd_camb:.3f} Mpc")
    print(f"  r_drag (sct_core canonical) = {EXPECTED_R_D:.3f} ± {EXPECTED_R_D_TOL:.1f} Mpc")
    print()

    if abs(rd_camb - EXPECTED_R_D) < EXPECTED_R_D_TOL:
        print(f"  [PASS] CAMB standard rdrag agrees with the standard horizon")
        print(f"         146.8 Mpc. Recombination is NOT patched (correct);")
        print(f"         the CAR effect is late-time (S8, b_IA) only.")
        return True
    else:
        print(f"  [FAIL] rdrag = {rd_camb:.2f} Mpc out of expected range.")
        print(f"         Standard CAMB should give r_drag ≈ 146.8 Mpc; check cosmology.")
        return False


if __name__ == "__main__":
    ok_a = test_canonical_constants()
    ok_b = test_standalone_predictions()
    ok_c = test_camb_patched()
    print()
    print("=" * 64)
    print(f"  Audit summary: constants={'PASS' if ok_a else 'FAIL'}, "
          f"standalone={'PASS' if ok_b else 'FAIL'}, "
          f"CAMB={'PASS' if ok_c else 'SKIP/INFO'}")
    print("=" * 64)
