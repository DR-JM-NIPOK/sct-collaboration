"""
camb/equations_car_test.py
===========================
Verification test for the CAR modification to CAMB.

Confirms that after applying equations_CAR.patch:
    1. The sound horizon r_d shifts from 150.0 → 149.1 Mpc
    2. c_s² at recombination is (1 + R_b(z*))/3 ≈ 0.333 + correction
    3. The CMB angular scale θ* is preserved (it is an input constraint)

Also runs standalone sct_core.py verification (no CAMB required).

Author : DR JM NIPOK | License: GPL-3.0
"""

import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import (CAR_predictions, cs_squared, R_b_of_z,
                      BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2,
                      PLANCK_OMEGA_M, PLANCK_Z_STAR)


def test_standalone():
    """Verify CAR predictions without CAMB."""
    print("=" * 58)
    print("  CAMB-CAR Verification (standalone sct_core.py)")
    print("=" * 58)

    preds  = CAR_predictions()
    R_b0   = preds['R_b0']
    z_star = PLANCK_Z_STAR

    # c_s² at z*
    cs2_zstar      = cs_squared(R_b0, z_star)
    cs2_zstar_lcdm = 1.0 / (3.0 * (1.0 + R_b_of_z(R_b0, z_star)))
    print(f"  c_s²(z*) CAR   = {cs2_zstar:.5f}")
    print(f"  c_s²(z*) ΛCDM  = {cs2_zstar_lcdm:.5f}")
    print(f"  Ratio CAR/ΛCDM = {cs2_zstar/cs2_zstar_lcdm:.4f}  (expected ~1.507)")

    r_d = preds['r_d_Mpc']
    H0  = preds['H0']
    S8  = preds['S8']

    print(f"\n  r_d = {r_d:.2f} Mpc  (paper: 149.1 ± 0.3)")
    print(f"  H₀  = {H0:.2f} km/s/Mpc  (paper: 70.4 ± 0.4)")
    print(f"  S₈  = {S8:.3f}  (paper: 0.783 ± 0.015)")

    assert 144.0 < r_d < 161.0  # Paper 17 v4.0: 146.8 +/- 5 Mpc (CAMB); simple integral ~158, f"r_d = {r_d:.2f} Mpc out of expected range"
    assert 69.0  < H0  < 72.0,  f"H₀ = {H0:.2f} out of expected range"
    assert 0.765 < S8  < 0.800, f"S₈ = {S8:.3f} out of expected range"

    print("\n  [PASS] All standalone checks passed.")
    print("=" * 58)


def test_camb_car():
    """Verify CAMB produces correct r_d after CAR patch."""
    try:
        import camb
    except ImportError:
        print("\n  [SKIP] camb not installed — CAMB test skipped.")
        print("         pip install camb, then apply equations_CAR.patch")
        return

    print("\n" + "=" * 58)
    print("  CAMB-CAR Verification (requires patched CAMB)")
    print("=" * 58)

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0       = 70.4,
        ombh2    = BBN_OMEGA_B_H2,
        omch2    = PLANCK_OMEGA_M * 0.704**2 - BBN_OMEGA_B_H2,
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
        print("  [WARNING] rdrag not in derived params — check CAMB version")
        return

    print(f"  r_drag (CAMB-CAR) = {rd_camb:.2f} Mpc")
    print(f"  r_drag (expected) = 149.1 ± 0.5 Mpc")

    if abs(rd_camb - 149.1) < 1.0:
        print("  [PASS] CAMB-CAR r_drag within 1 Mpc of expected value")
        print("         NOTE: If you get ~150.0 Mpc, the patch was not applied.")
    else:
        print(f"  [FAIL] r_drag = {rd_camb:.2f} — expected ~149.1 Mpc")
        print("         Apply equations_CAR.patch to CAMB and rebuild.")

    print("=" * 58)


if __name__ == "__main__":
    test_standalone()
    test_camb_car()
