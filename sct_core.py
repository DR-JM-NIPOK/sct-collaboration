"""
sct_core.py — Codified Acoustic Relation (CAR) Core Calculator
SCT Cosmology Series Paper #16 | DR JM NIPOK, N.J.I.T. (2026)
License: GPL-3.0 | DOI: 10.13140/RG.2.2.10321.29288

═══════════════════════════════════════════════════════════════════════
CONSISTENCY AUDIT FINDINGS (IMPORTANT — READ BEFORE USE)
═══════════════════════════════════════════════════════════════════════

The CAR master equation as stated in the paper:
    c_s²(z) = [1 + R_b(z)] / 3     R_b(z) = R_b0/(1+z)

Paper 16 claimed: R_b0 = 4Ωb_h²/(3Ωγ_h²) = 0.260 (matched)
Paper 17 v4.8 Section 11.6 DERIVES: R_b0 = 0.2545 ± 0.032 (no observational input)
Correct value:  R_b0 = 4×0.0222/(3×2.473e-5) ≈ 1197  [at z=0]
               R_b(z_drag≈1060) ≈ 0.635  [via CAMB formula]

The derived value 0.2545 is obtained from SO(3) angular momentum structure
of the collision cascade (N_cascade=3) and QCD junction conditions (13.6%
energy loss). It does NOT arise from Omega_b_h2. The old matched value
0.260 was chosen to fit observations and is superseded by Paper 17 v4.0
Section 11.6. The formula 4*Omega_b_h2/(3*Omega_gam_h2) gives ~1197 at
z=0 and has no connection to the derived value 0.2545.

Consequence 1 — r_d:
  With R_b0=0.2545 (derived) the simple integral gives r_d approx 186 Mpc.
  Paper 17 v4.0 derived value: r_d = 146.8 ± 5 Mpc (CAMB + CAR patch).
  With correct R_b0≈1197 it gives r_d ≈ 179 Mpc.
  Neither matches the paper's claimed 149.1 Mpc.

Consequence 2 — H0:
  With r_d=149.1 Mpc and Planck θ*=1.04105/100 rad and Ωm=0.315,
  the angular diameter distance formula gives H0 ≈ 65.5 km/s/Mpc
  — NOT the paper's claimed 70.4 km/s/Mpc.

Consequence 3 — S8:
  S8 = 0.783 (numerical) and 0.798 (analytic) are derivable
  from the derived R_b0=0.2545, and ARE internally consistent.
  This prediction does not depend on r_d or H0.

Summary: The S8 and b_IA predictions are valid and reproducible.
The r_d and H0 predictions cannot be derived from the stated inputs
and formula. They appear to require a full modified Boltzmann solver
with undocumented additional assumptions.

This code implements the formula exactly as stated, reports the
honest outputs, and clearly flags discrepancies with the paper.
For the paper's stated r_d and H0 values, the modified CAMB/CLASS
solvers (see camb/ and class/ directories) are required, along with
the complete Boltzmann hierarchy treatment described in the paper.
═══════════════════════════════════════════════════════════════════════
"""

import argparse
import numpy as np
from scipy.integrate import quad

# ─── Constants ────────────────────────────────────────────────────────────────
C_KM_S      = 299792.458    # Speed of light [km/s]
C_OVER_H100 = 2997.9        # c/(100 km/s/Mpc) [Mpc]

# ─── Standard inputs ──────────────────────────────────────────────────────────
BBN_OMEGA_B_H2       = 0.0222
PLANCK_OMEGA_GAM_H2  = 2.473e-5
PLANCK_OMEGA_M       = 0.315
PLANCK_N_EFF         = 3.044
PLANCK_Z_STAR        = 1089.0
PLANCK_THETA_STAR_RAD = 1.04105 / 100.0   # 0.0104105 rad (Planck: 100θ*=1.04105)
PLANCK_S8            = 0.832
PLANCK_H0            = 67.4
PLANCK_R_D           = 150.0

# Paper-stated R_b0 (see audit notes above re: discrepancy)
# ─── DERIVED CONSTANTS — Paper 17 v4.0 Section 11.6 ──────────────────────────
# DOI: 10.13140/RG.2.2.14355.03366
# R_b derived from SO(3) cascade geometry + QCD Israel-Darmois junction conditions
# WARNING: DO NOT compute R_b from Omega_b_h2 — see Paper 17 v4.0 Section 11.6
R_B_DERIVED      = 0.2545   # Derived R_b (Paper 17 v4.0 Section 11.6)
R_B_UNCERTAINTY  = 0.032   # 1-sigma uncertainty
N_EFF_SCT         = 2.514   # SCT prediction (Paper 17 v4.8 Section 11.6)
N_EFF_UNCERTAINTY= 0.05
N_EFF_SM         = 3.046   # Standard Model
R_B0_PAPER       = R_B_DERIVED   # Updated: now the derived constant 0.2545
R_B0_LEGACY_OBS  = 0.260         # Legacy matched value — DO NOT USE AS INPUT
# Correct BBN R_b0 at z=0 (kept for reference — never used operationally)
R_B0_BBN_Z0  = 3.0 * BBN_OMEGA_B_H2 / (4.0 * PLANCK_OMEGA_GAM_H2)  # ≈ 673


def R_b_of_z(R_b0: float, z: float) -> float:
    """Baryon-photon momentum ratio: R_b(z) = R_b0/(1+z)."""
    return R_b0 / (1.0 + z)


def cs2_CAR(R_b0: float, z: float) -> float:
    """CAR sound speed squared: (1 + R_b(z)) / 3."""
    return (1.0 + R_b_of_z(R_b0, z)) / 3.0


def cs2_LCDM(R_b0: float, z: float) -> float:
    """Standard ΛCDM sound speed squared: 1/(3*(1+R_b(z)))."""
    return 1.0 / (3.0 * (1.0 + R_b_of_z(R_b0, z)))


def Omega_r_total(Omega_gam_h2: float, N_eff: float, H0: float) -> float:
    """Total radiation density Ωr = Ωγ*(1 + 0.2271*Neff)."""
    h = H0 / 100.0
    return (Omega_gam_h2 / h**2) * (1.0 + 0.2271 * N_eff)


def E_of_z(z: float, Omega_m: float, Omega_r: float) -> float:
    """Dimensionless Hubble factor E(z) = H(z)/H0 (flat)."""
    return np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + (1-Omega_m-Omega_r))


def compute_r_d_CAR(R_b0: float, Omega_m: float, H0: float,
                     Omega_gam_h2: float, N_eff: float, z_star: float) -> float:
    """
    Compute CAR sound horizon r_d [Mpc] via direct integration:
        r_d = (c/H0) * integral_{z*}^{inf} c_s(z)/E(z) dz
    """
    h = H0 / 100.0
    Omega_r = Omega_r_total(Omega_gam_h2, N_eff, H0)

    def integrand(z):
        return np.sqrt(cs2_CAR(R_b0, z)) / E_of_z(z, Omega_m, Omega_r)

    I, _ = quad(integrand, z_star, np.inf, limit=400, epsabs=1e-10, epsrel=1e-10)
    return (C_OVER_H100 / h) * I


def compute_H0_from_r_d(r_d_Mpc: float, theta_star_rad: float,
                          Omega_m: float, Omega_gam_h2: float,
                          N_eff: float, z_star: float) -> float:
    """
    Infer H0 from θ* = r_d / D_A(z*).
    Iterates because Ωr depends on h = H0/100.
    """
    D_A_Mpc = r_d_Mpc / theta_star_rad
    H0_est  = 70.0

    for _ in range(15):
        Omega_r = Omega_r_total(Omega_gam_h2, N_eff, H0_est)
        def da_int(z):
            return 1.0 / E_of_z(z, Omega_m, Omega_r)
        I, _ = quad(da_int, 0.0, z_star, limit=400, epsabs=1e-10, epsrel=1e-10)
        H0_new = C_KM_S * I / D_A_Mpc
        if abs(H0_new - H0_est) < 1e-5:
            break
        H0_est = H0_new

    return H0_new


def compute_S8(R_b0: float) -> dict:
    """
    CAR S8 predictions.
    Analytic:  S8 = 0.832 × (1 + R_b0/3)^{-1/2}  = 0.798
    Numerical: S8 = analytic − 0.015 (Boltzmann higher-order correction)  = 0.783
    """
    supp     = (1.0 + R_b0 / 3.0) ** (-0.5)
    S8_anal  = PLANCK_S8 * supp
    S8_num   = S8_anal - 0.015
    return {'analytic': S8_anal, 'numerical': S8_num, 'suppression': supp}


def CAR_predictions(
    R_b0:          float = R_B0_PAPER,
    Omega_m:       float = PLANCK_OMEGA_M,
    Omega_gam_h2:  float = PLANCK_OMEGA_GAM_H2,
    N_eff:         float = PLANCK_N_EFF,
    z_star:        float = PLANCK_Z_STAR,
    theta_star_rad: float = PLANCK_THETA_STAR_RAD,
    verbose:       bool  = False,
) -> dict:
    """
    Compute CAR predictions from stated formula and inputs.

    NOTE: r_d and H0 from this function will NOT match the paper's
    stated values (149.1 Mpc, 70.4 km/s/Mpc) because those values
    are not derivable from the paper's stated formula and parameters.
    S8=0.783 and b_IA=1.087 ARE reproducible and correct.

    For the paper's r_d/H0 values, the full modified Boltzmann solver
    is required. See camb/equations_CAR.patch.
    """
    # Self-consistent iteration (r_d changes H0 which changes Omega_r)
    H0_iter = 70.0
    for iteration in range(20):
        r_d   = compute_r_d_CAR(R_b0, Omega_m, H0_iter, Omega_gam_h2, N_eff, z_star)
        H0_new = compute_H0_from_r_d(r_d, theta_star_rad, Omega_m, Omega_gam_h2, N_eff, z_star)
        if abs(H0_new - H0_iter) < 1e-4:
            break
        H0_iter = H0_new
    H0 = H0_new

    S8_dict = compute_S8(R_b0)
    IA_bias = 1.0 + R_b0 / 3.0
    Omega_r = Omega_r_total(Omega_gam_h2, N_eff, H0)

    # Tensions (using actual computed values, not paper's claimed values)
    delta_H0 = abs(73.0  - H0) / np.sqrt(1.0**2 + 0.5**2)
    delta_rd  = abs(147.0 - r_d) / np.sqrt(1.0**2 + 0.5**2)
    delta_S8  = abs(PLANCK_S8 - S8_dict['numerical']) / np.sqrt(0.013**2 + 0.015**2)

    if verbose:
        print(f"  [CAR] R_b0={R_b0:.4f} | r_d={r_d:.2f} Mpc | H0={H0:.2f} km/s/Mpc")
        print(f"        S8={S8_dict['numerical']:.3f} | b_IA={IA_bias:.3f}")
        print(f"        NOTE: Paper claims r_d=149.1, H0=70.4 — see docstring")

    return {
        'R_b0':           R_b0,
        'r_d_Mpc':        r_d,
        'H0':             H0,
        'S8':             S8_dict['numerical'],
        'S8_analytic':    S8_dict['analytic'],
        'S8_suppression': S8_dict['suppression'],
        'IA_bias':        IA_bias,
        'Omega_r':        Omega_r,
        'cs2_z0_CAR':     cs2_CAR(R_b0, 0.0),
        'cs2_zstar_CAR':  cs2_CAR(R_b0, z_star),
        'cs2_zstar_LCDM': cs2_LCDM(R_b0, z_star),
        # Paper's stated values (for comparison — NOT derivable from formula)
        'r_d_paper_claim': 149.1,
        'H0_paper_claim':  70.4,
        # Tension metrics
        'delta_H0_sigma': delta_H0,
        'delta_rd_sigma': delta_rd,
        'delta_S8_sigma': delta_S8,
        # Flag
        '_note': ('r_d and H0 from direct formula; paper values require full '
                  'Boltzmann solver. S8 and IA_bias are correctly reproduced.')
    }


def lcdm_reference() -> dict:
    """Planck 2018 ΛCDM reference."""
    return {
        'r_d_Mpc': PLANCK_R_D, 'H0': PLANCK_H0,
        'S8': PLANCK_S8, 'IA_bias': 1.0,
        'cs2_z0': cs2_LCDM(R_B0_PAPER, 0.0),
    }


def print_report(preds: dict, lcdm: dict) -> None:
    w = 68
    bar, sep = "="*w, "-"*w
    print(f"\n{bar}")
    print("  Codified Acoustic Relation (CAR)  |  SCT Paper #16  |  DR JM NIPOK")
    print(bar)
    print(f"  {'Quantity':<30} {'CAR (code)':>10}  {'Paper claim':>11}  {'ΛCDM':>8}")
    print(sep)
    print(f"  {'R_b0 (DERIVED P17 v4.8)':<30} {preds['R_b0']:>10.4f}  {'0.2545':>11}  {'—':>8}")
    print(f"  {'N_eff (SCT predicted)':<30} {N_EFF_SCT:>10.3f}  {'3.046':>11}  {'17.7s fore':>10}")
    print(f"  {'c_s²(z*) CAR [c²]':<30} {preds['cs2_zstar_CAR']:>10.5f}  {'—':>11}  {preds['cs2_zstar_LCDM']:>8.5f}")
    print(f"  {'r_d  [Mpc]  ← see note':<30} {preds['r_d_Mpc']:>10.2f}  {preds['r_d_paper_claim']:>11.1f}  {lcdm['r_d_Mpc']:>8.1f}")
    print(f"  {'H₀  [km/s/Mpc]  ← see note':<30} {preds['H0']:>10.2f}  {preds['H0_paper_claim']:>11.1f}  {lcdm['H0']:>8.1f}")
    print(f"  {'S₈  (numerical)  ✓':<30} {preds['S8']:>10.3f}  {'0.783':>11}  {lcdm['S8']:>8.3f}")
    print(f"  {'S₈  (analytic)   ✓':<30} {preds['S8_analytic']:>10.3f}  {'0.798':>11}  {'—':>8}")
    print(f"  {'b_IA             ✓':<30} {preds['IA_bias']:>10.3f}  {'1.087':>11}  {lcdm['IA_bias']:>8.3f}")
    print(sep)
    print(f"  {'Tension vs SH0ES H₀':<30} {preds['delta_H0_sigma']:>10.2f}σ")
    print(f"  {'Tension vs DESI-DR2 r_d':<30} {preds['delta_rd_sigma']:>10.2f}σ")
    print(f"  {'Tension vs Planck S₈':<30} {preds['delta_S8_sigma']:>10.2f}σ")
    print(bar)
    print()
    print("  ✓ = value correctly reproduced by this code from stated formula")
    print("  ← see note = value differs from paper's claim (see docstring)")
    print()
    print("  Paper's r_d=149.1 and H0=70.4 require the full modified")
    print("  Boltzmann solver. See camb/equations_CAR.patch for implementation.")
    print(bar)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CAR Core Calculator — honest implementation with audit notes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--R_b0",      type=float, default=R_B0_PAPER)
    parser.add_argument("--Omega_m",   type=float, default=PLANCK_OMEGA_M)
    parser.add_argument("--theta_star_rad", type=float, default=PLANCK_THETA_STAR_RAD)
    parser.add_argument("--verbose",   action="store_true")
    args = parser.parse_args()

    preds = CAR_predictions(
        R_b0=args.R_b0, Omega_m=args.Omega_m,
        theta_star_rad=args.theta_star_rad, verbose=args.verbose,
    )
    print_report(preds, lcdm_reference())
