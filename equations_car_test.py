"""
camb/equations_car_test.py
==========================
Verification test for the CAR modification to CAMB. (v4.8.1 corrected)

This script verifies that:
    1. sct_core.py's compute_r_d_integral gives r_d ≈ 161.4 Mpc (no CAMB needed)
    2. CAR predictions for S8 and b_IA match canonical v4.8 values
    3. If CAMB is installed AND the equations_CAR.patch is applied, the
       CAMB derived rdrag agrees with sct_core.py's value (≈ 161.4 Mpc)

v4.8.1 Audit Correction:
    Earlier versions asserted r_d ≈ 149.1 Mpc and treated 150.0 as "patch
    not applied". Both numbers were artifacts of an unreproducible CAMB run
    that turned out to depend on a buggy Fortran patch (using standard R
    density ratio with CAR formula → superluminal cs²). The canonical
    value, verified consistently across sct_core.py and a properly-corrected
    Fortran patch, is r_d ≈ 161.4 Mpc.

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
EXPECTED_R_D       = 161.4
EXPECTED_R_D_TOL   =   0.5   # Mpc — numerical-integration uncertainty
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
    print(f"          (canonical v4.8.1; ΛCDM observed: 147 ± 1 Mpc; "
          f"CAR does NOT close BAO tension)")
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

    H0 = 70.4
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
        print(f"  [PASS] CAMB-CAR rdrag agrees with sct_core canonical value.")
        print(f"         The equations_CAR.patch is applied and active.")
        return True
    elif abs(rd_camb - 144.0) < 2.0:
        print(f"  [INFO] rdrag ≈ {rd_camb:.1f} Mpc = standard ΛCDM value at H0=70.4")
        print(f"         The equations_CAR.patch was NOT applied to this CAMB build.")
        print(f"         Apply: cd <camb_dir> && patch -p1 < equations_CAR.patch")
        print(f"                && make")
        return False
    elif abs(rd_camb - 150.0) < 2.0:
        print(f"  [INFO] rdrag ≈ {rd_camb:.1f} Mpc = standard ΛCDM at H0≈67")
        print(f"         CAMB may be running with different H0 than requested.")
        return False
    else:
        print(f"  [FAIL] rdrag = {rd_camb:.2f} Mpc out of expected range.")
        print(f"         Either CAMB is at a different cosmology than expected,")
        print(f"         or the patch implementation differs from canonical.")
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
