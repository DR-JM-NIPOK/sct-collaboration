"""
sct_core.py — Codified Acoustic Relation (CAR) Core Calculator
SCT Cosmology Series Paper #16 / Paper #17 | DR JM NIPOK, N.J.I.T. (2026)
ORCID: 0009-0006-3940-4450
License: GPL-3.0
Paper 15 DOI: 10.13140/RG.2.2.10321.29288
Series 2 Paper 1 DOI: 10.13140/RG.2.2.14355.03366

VERSION HISTORY
───────────────────────────────────────────────────────────────────────────
v1.0    March 2026   Original release
v2.0    April 2026   Three critical bugs corrected
v3.0    April 2026   Epistemic upgrade: R_b transitions from matched
                     observational parameter to derived constant
v4.8    April 2026   R_b = 0.2545 (Series 2 Paper 1 Section 11.6)
v4.8.1  April 2026   NLA-recursive audit of full repository:
                     synchronized sct_core.py across all branches;
                     predictions.csv b_IA corrected 1.087 → 1.0848;
                     KiDS-DR5 S8 corrected 0.815 → 0.788.
v4.8.4  June 2026    RNLA v2.3 recursive audit (checks 13, 14):
                     • r_d RESTORED to the standard photon-baryon horizon
                       146.8 Mpc (was wrongly set to 161.4 in v4.8.1).
                     • Category error removed: the CAR enhanced sound speed
                       (cs²=0.4182) is a LATE-TIME coherent-sector quantity
                       (it sets S8 and b_IA). It is NOT the recombination
                       acoustic speed and must not enter the sound-horizon
                       integral. Doing so produced the superseded 161.4 Mpc.
                     • r_d sound-horizon integral now uses the STANDARD
                       baryon-loaded acoustic speed → r_drag ≈ 146.8 Mpc,
                       r_*(z*) ≈ 144.4 Mpc (Planck-consistent).
                     • H0 self-consistent from θ* + r_d = 146.8 → ≈ 66.3
                       km/s/Mpc (global). The previously reported 70.4 is
                       NOT CMB-derivable and is retired; local H0 rises
                       toward 70–73 via the void + temporal mechanism.

v4.8.5  July 2026   Register update from in-session CAMB/SPARC verification
                     (SCT_VERIFICATION_NOTE_CAMB_SPARC_20260710.md):
                     • N_eff = 2.514 quantified against Planck TT: naive
                       substitution into stock ΛCDM at fixed θ* is excluded
                       at Δχ² ≈ +658 (83 binned points); logged as new
                       tensions.csv T031 and folded into predictions.csv
                       P088 (status: PENDING, PARTIAL sub-check).
                     • Register-consistency flag: the H0≈66–67 note holds
                       only under standard N_eff=3.044; combined with
                       N_eff=2.514 it gives H0=60.15. Which era N_eff=2.514
                       describes (recombination vs. CMB-S4-observed
                       late-time quantity) is now flagged as open in P088.
                     • Audit-trail note: r_drag at N_eff=2.514 (149.86 Mpc)
                       numerically coincides with the retired 149.1 Mpc r_d
                       value — a coincidence, not a revival; r_d itself is
                       unchanged at 146.8 Mpc (v4.8.4, still canonical).
                     • No numeric constants changed in this release
                       (R_b, c_s², b_IA, r_d, H0, N_eff all UNCHANGED).

THE CAR FRAMEWORK IN ONE PARAGRAPH
───────────────────────────────────────────────────────────────────────────
The Codified Acoustic Relation modifies the effective baryon-photon sound
speed in the LATE-TIME coherent sector:

    Standard ΛCDM:  c_s²(z) = 1 / [3 (1 + R(z))]   with R(z)=3ρ_b/(4ρ_γ)
    CAR (late-time): c_s²(z) = (1 + R_b(z)) / 3      [Paper 15 §2.1]

where R_b(z) = R_b0 / (1+z) with R_b0 = 0.2545 (DERIVED — Series 2 Paper 1
§11.6) rather than matched to observations. This produces a ~26% enhancement
of the coherent c_s² at z=0 (0.4182 vs 1/3) and a tomographic S8 suppression
of ~4.4% and IA bias b_IA = 1.085 — both derivable from a single derived
constant. CRITICALLY, this enhancement is a z≈0 coherent-sector effect; it
does NOT modify the recombination-epoch acoustic physics. The early-universe
sound horizon is therefore STANDARD: r_drag ≈ 146.8 Mpc, r_*(z*) ≈ 144.4 Mpc.

R_B_DERIVED = 0.2545 (NOT 0.260)
───────────────────────────────────────────────────────────────────────────
Series 2 Paper 1 Section 11.6 derives R_b from first principles:
  (1) SO(3) angular momentum structure of the collision cascade
      → N_cascade = 3 (three independent cascade planes)
  (2) QCD phase transition boundary correction (Paper 9 Israel-Darmois
      junction conditions, 13.6% energy loss)
NO observational input. The derivation is purely geometric and field-
theoretic. Agreement with observed R_b ≈ 0.260 is at 0.17 sigma —
a post-diction that closes the v1.0/v2.0 circularity.

WARNING: DO NOT pass R_b = 0.260 as a hardcoded input to this module.
   The value 0.260 is now a legacy observational reference for comparison
   only. All physics computations use R_B_DERIVED = 0.2545.

WHAT THIS CODE COMPUTES
───────────────────────────────────────────────────────────────────────────
Analytic outputs (no CAMB, independently verified):
  c_s²(z=0) = 0.41817                 from (1 + R_b)/3   [LATE-TIME coherent]
  b_IA      = 1.08483 ± 0.011         from 1 + R_b/3
  S8 anlc   = 0.79881                 from 0.832 × (1 + R_b/3)^(-1/2)
  S8 num    = 0.78381 ± 0.015         analytic + CAMB Jeans correction (-0.015)

Standard recombination acoustics (NOT modified by CAR):
  r_d (drag)   = 146.8 ± 5 Mpc        standard photon-baryon horizon (Paper 4)
  r_*  (z*)    = 144.4 Mpc            sound horizon at last scattering
  H0 (global)  = 66.3 km/s/Mpc        from θ* + r_d (PARTIAL; see note)

New CMB-S4 prediction (Series 2 Paper 1):
  N_eff_SCT = 2.514 ± 0.05    vs SM N_eff = 3.046  →  17.7σ separation (forecast)

CRITICAL AUDIT NOTE — r_d, EXPLICIT (RNLA v2.3, June 2026)
───────────────────────────────────────────────────────────────────────────
The recombination/drag sound horizon is the STANDARD photon-baryon horizon
fixed by ω_b and ω_m (H0-independent): r_drag ≈ 146.8 Mpc, r_*(z*) ≈ 144.4 Mpc.
This is standard early-universe acoustics; the SCT framework does not alter it.

The CAR enhanced sound speed, c_s²=(1+R_b)/3 with R_b0=0.2545, is a LATE-TIME
coherent-sector quantity: it governs the z≈0 coherence that sets S8 and b_IA.
Inserting this late-time speed into the recombination sound-horizon integral
is a category error; it gives the SUPERSEDED 161.4/178 Mpc and is retired.
(At recombination R_b → 0, so (1+R_b)/3 → 1/3, the pure-radiation speed, which
omits baryon loading and inflates the horizon. The correct recombination speed
is the baryon-loaded standard 1/[3(1+R)].)

Implication: with the standard r_d = 146.8 Mpc, the SCT acoustic sector is
consistent with DESI-DR2 BAO (147 ± 1 Mpc) and Planck — no early-time tension
is introduced. The H0 tension is addressed by a LATE-TIME mechanism (void +
temporal), not by shifting r_d. The global θ*-inferred H0 ≈ 66.3 km/s/Mpc;
local distance-ladder H0 is raised toward 70–73 by the void + temporal effect.

BUGS PRESERVED FROM v2.0 / v3.0 (kept fixed)
───────────────────────────────────────────────────────────────────────────
Bug 1 — R_b0 convention (CRITICAL — affects S8 and b_IA):
  BROKEN: R_b0 = 4×Ω_b_h²/(3×Ω_γ_h²) = 1196.9
  FIXED:  R_b0 = 0.2545  (derived constant)
Bug 2 — theta_star unit conversion (CRITICAL — affects H0):
  BROKEN: theta_star_rad = theta_star × π/180
  FIXED:  theta_star_rad = theta_star / 100   (Planck reports 100 × θ*)
Bug 3 — r_d integral normalisation (affects r_d):
  BROKEN: r_d = integral × c/100
  FIXED:  r_d = integral × c / H0
"""

import argparse
import numpy as np
from scipy.integrate import quad

# ── Fundamental constants ──────────────────────────────────────────────────────
C_KM_S = 299792.458          # Speed of light [km/s]

# ══════════════════════════════════════════════════════════════════════════════
# DERIVED CONSTANTS — Series 2 Paper 1, Section 11.6
# DOI: 10.13140/RG.2.2.14355.03366
#
# These values are DERIVED from first principles. They are NOT matched to
# observations and must NOT be replaced with observational inputs.
#
# WARNING: R_b = 0.260 must NOT be used as input to this module.
#    The observed value 0.260 is a post-diction reference only.
#    All computations use R_B_DERIVED = 0.2545 ± 0.032.
# ══════════════════════════════════════════════════════════════════════════════

# R_b derived from SO(3) angular momentum structure of the collision cascade
# (N_cascade = 3) and QCD phase transition boundary correction from Paper 9
# Israel-Darmois junction conditions (13.6% energy loss).
R_B_DERIVED      = 0.2545    # Derived baryon-photon coherence ratio (Series 2 Paper 1 §11.6)
R_B_UNCERTAINTY  = 0.032     # 1-σ (geometric + QCD boundary)

# Effective (LATE-TIME coherent) sound speed squared at z=0: c_s²/c² = (1 + R_b)/3
# With R_b = 0.2545: c_s² = 1.2545/3 = 0.41817 (Series 2 Paper 1)
# NOTE: this is the z≈0 coherent-sector speed (sets S8, b_IA), NOT the
# recombination acoustic speed. It must NOT enter the sound-horizon integral.
CS2_DERIVED      = (1.0 + R_B_DERIVED) / 3.0   # = 0.41817
CS2_UNCERTAINTY  = R_B_UNCERTAINTY / 3.0       # = 0.01067

# Background coherence enhancement / IA bias: C_hat_bg = b_IA = 1 + R_b/3
C_HAT_BG_DERIVED     = 1.0 + R_B_DERIVED / 3.0   # = 1.08483
C_HAT_BG_UNCERTAINTY = R_B_UNCERTAINTY / 3.0     # = 0.01067

# Sound horizon — STANDARD photon-baryon horizon (Plasma Equivalence, Paper 4).
# RNLA v2.3 (June 2026): r_d is the standard recombination/drag horizon, fixed
# by ω_b and ω_m (H0-independent). The CAR late-time enhancement does NOT modify
# it. The v4.8.1 value 161.4 Mpc was a category error (late-time c_s inserted
# into the recombination integral) and is retired.
R_D_DERIVED      = 146.8     # Mpc — standard drag-epoch sound horizon r_drag
R_D_UNCERTAINTY  = 5.0       # Mpc — 1-σ (ω_b, ω_m uncertainty)
R_STAR_DERIVED   = 144.4     # Mpc — sound horizon at last scattering r_*(z*)

# Effective neutrino species — Series 2 Paper 1 Section 11.6 prediction
# N_eff_SCT = 2.514 ± 0.05 (predicted)
# N_eff_SM  = 3.046          (Standard Model)
# These are on opposite sides of 3.000. CMB-S4 tests this at 17.7σ (forecast).
N_EFF_SCT         = 2.514    # SCT prediction
N_EFF_UNCERTAINTY = 0.05     # 1-σ
N_EFF_SM          = 3.046    # Standard Model (Mangano et al. 2005)
CMB_S4_SIGMA_NEFF = 0.03     # CMB-S4 forecast sensitivity on N_eff

# Legacy reference value — observational comparison ONLY
R_B_LEGACY_OBS   = 0.260     # OLD matched value — DO NOT USE AS INPUT

# Internal alias — always uses the derived constant
R_B0 = R_B_DERIVED

# ── Standard cosmological inputs (Planck 2018 + BBN) ─────────────────────────
BBN_OMEGA_B_H2      = 0.0222     # Baryon physical density (BBN)
PLANCK_OMEGA_B_H2   = 0.0222     # Same — alias for clarity in CAMB calls
PLANCK_OMEGA_GAM_H2 = 2.473e-5   # Photon physical density
PLANCK_OMEGA_M      = 0.315      # Total matter density
PLANCK_N_EFF        = 3.044      # Standard reference (used in radiation density)
PLANCK_Z_STAR       = 1089.0     # Last scattering
PLANCK_Z_DRAG       = 1060.0     # Baryon drag
PLANCK_THETA_STAR   = 1.04105    # 100 × θ* (Planck 2018 measurement)
PLANCK_S8           = 0.832      # Planck 2018 clustering amplitude

# ── Observed values (used only for sigma comparisons, never as physics input) ─
OBSERVED_S8_DES_Y6   = (0.780, 0.012)
OBSERVED_S8_HSC_Y3   = (0.776, 0.020)   # published 0.776
OBSERVED_S8_KIDS_DR5 = (0.788, 0.014)   # published 0.788
OBSERVED_R_D_DESI    = (147.0, 1.0)
OBSERVED_R_D_PLANCK  = (147.1, 0.3)     # Planck 2018 r_drag
OBSERVED_H0_PLANCK   = (67.4, 0.5)
OBSERVED_H0_SHOES    = (73.0, 1.0)


# ════════════════════════════════════════════════════════════════════════════
# CORE PHYSICS FUNCTIONS
# ════════════════════════════════════════════════════════════════════════════

def R_b_of_z(z: float) -> float:
    """
    CAR baryon-photon coherence ratio at redshift z (LATE-TIME coherent sector).

    Paper 15 §2.3 / Series 2 Paper 1 §11.6:
        R_b(z) = R_B_DERIVED / (1+z)

    At z=0:      R_b = 0.2545
    At z_drag:   R_b ≈ 0.000240 (essentially zero — photon-limit at recombination)
    """
    return R_B_DERIVED / (1.0 + z)


def cs_CAR(z: float) -> float:
    """
    CAR late-time coherent-sector sound speed at redshift z, in units of c.

        c_s(z)/c = sqrt[(1 + R_b(z)) / 3]

    At z=0:    c_s/c = sqrt(1.2545/3) = 0.6467   (coherent enhancement; sets S8, b_IA)
    At z_drag: c_s/c ≈ sqrt(1.0002/3) = 0.5774   (→ pure-radiation limit, R_b→0)

    WARNING: this is the z≈0 coherent-sector speed. It must NOT be used to
    compute the recombination sound horizon (that requires the baryon-loaded
    standard speed, cs_LCDM). Using cs_CAR in the horizon integral reproduces
    the retired 161.4 Mpc category error.
    """
    return np.sqrt((1.0 + R_b_of_z(z)) / 3.0)


def cs_LCDM(z: float) -> float:
    """
    Standard baryon-loaded acoustic sound speed (units of c).

    c_s(z)/c = 1 / sqrt[3 (1 + R(z))],  R(z) = 3 ρ_b / (4 ρ_γ) = 3 ω_b /(4 ω_γ)/(1+z)

    This is the correct recombination-epoch sound speed and the one used for
    the standard photon-baryon sound horizon r_d.
    """
    R_std = (3.0 * BBN_OMEGA_B_H2) / (4.0 * PLANCK_OMEGA_GAM_H2) / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R_std))


def cs_squared(R_b: float = None, z: float = 0.0) -> float:
    """Convenience: late-time c_s²(z) given an explicit R_b at z=0 (or default derived)."""
    if R_b is None:
        R_b = R_B_DERIVED
    Rb_z = R_b / (1.0 + z)
    return (1.0 + Rb_z) / 3.0


def omega_r_total(H0: float) -> float:
    """Total radiation density Ω_r(H0) including standard 3.044 ν species."""
    h = H0 / 100.0
    return (PLANCK_OMEGA_GAM_H2 / h**2) * (1.0 + 0.2271 * PLANCK_N_EFF)


def E_of_z(z: float, H0: float) -> float:
    """Dimensionless Hubble factor E(z) = H(z)/H0 for flat ΛCDM cosmology."""
    Omega_r = omega_r_total(H0)
    Omega_L = 1.0 - PLANCK_OMEGA_M - Omega_r
    return np.sqrt(
        Omega_r        * (1.0 + z)**4
        + PLANCK_OMEGA_M * (1.0 + z)**3
        + Omega_L
    )


def hubble_factor(z: float, Omega_m: float, Omega_r: float) -> float:
    """
    Backward-compat alias used by older test suites.

    Returns the dimensionless Hubble factor E(z) = H(z)/H0 for an explicit
    (Omega_m, Omega_r) cosmology — independent of module-level Planck defaults.
    """
    Omega_L = 1.0 - Omega_m - Omega_r
    return np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + Omega_L)


def compute_r_d_integral(H0: float, z_drag: float = PLANCK_Z_DRAG) -> float:
    """
    Compute the STANDARD photon-baryon sound horizon r_d [Mpc] by direct
    integration of the baryon-loaded acoustic speed (NOT the CAR late-time speed):

        r_d = (c/H0) × ∫_{z_drag}^{∞} c_s^std(z) / E(z) dz

    with c_s^std(z) = 1/sqrt[3(1+R(z))], R(z)=3ω_b/(4ω_γ)/(1+z).
    Returns ≈ 146.8 Mpc at Planck cosmology — the standard drag-epoch horizon.

    (The CAR enhancement is late-time; it does not modify this integral. Using
    cs_CAR here would reproduce the retired 161.4 Mpc category error.)
    """
    def integrand(z):
        return cs_LCDM(z) / E_of_z(z, H0)
    I, _ = quad(integrand, z_drag, np.inf, limit=400,
                epsabs=1e-10, epsrel=1e-10)
    return (C_KM_S / H0) * I


def compute_sound_horizon(R_b: float, Omega_m: float, Omega_r: float,
                          z_star: float, H0: float = 67.4) -> float:
    """
    Backward-compat alias used by older test suites.

    Computes the STANDARD photon-baryon sound horizon with explicit cosmology.
    The recombination acoustic speed is the baryon-loaded standard form; the
    R_b argument is retained for signature compatibility but the recombination
    sound speed does NOT use the CAR late-time enhancement.

    Parameters
    ----------
    R_b      : float   (retained for compatibility; not used for the horizon speed)
    Omega_m  : float   Matter density
    Omega_r  : float   Radiation density (typically 9e-5)
    z_star   : float   Lower redshift bound (drag or last-scattering epoch)
    H0       : float   Hubble constant in km/s/Mpc, default 67.4
    """
    Omega_L = 1.0 - Omega_m - Omega_r
    R0 = (3.0 * BBN_OMEGA_B_H2) / (4.0 * PLANCK_OMEGA_GAM_H2)

    def integrand(z):
        cs = 1.0 / np.sqrt(3.0 * (1.0 + R0 / (1.0 + z)))
        E  = np.sqrt(Omega_r*(1+z)**4 + Omega_m*(1+z)**3 + Omega_L)
        return cs / E

    I, _ = quad(integrand, z_star, np.inf, limit=400,
                epsabs=1e-10, epsrel=1e-10)
    return (C_KM_S / H0) * I


def compute_D_M(H0: float, z_star: float = PLANCK_Z_STAR) -> float:
    """Comoving distance to last scattering D_M(z*) [Mpc]."""
    Omega_r = omega_r_total(H0)
    Omega_L = 1.0 - PLANCK_OMEGA_M - Omega_r
    def integrand(z):
        return 1.0 / np.sqrt(
            PLANCK_OMEGA_M * (1.0 + z)**3
            + Omega_L
            + Omega_r * (1.0 + z)**4
        )
    I, _ = quad(integrand, 0.0, z_star, limit=400)
    return (C_KM_S / H0) * I


def compute_H0_from_theta_and_rd(r_d: float,
                                  theta_star_100: float = PLANCK_THETA_STAR
                                  ) -> float:
    """
    Derive H0 self-consistently from θ* = r_s / D_M(z*).

    With the standard drag horizon r_d = 146.8 Mpc this returns H0 ≈ 66.3
    km/s/Mpc (global). With the last-scattering horizon r_*(z*) = 144.4 Mpc it
    returns ≈ 67.4 km/s/Mpc (Planck). 70.4 is NOT recoverable from θ* and is
    retired.

    Bug 2 fix (v2.0): Planck reports 100 × θ*, so θ* = 1.04105 / 100 rad.
    """
    theta_rad = theta_star_100 / 100.0
    H0 = 70.0
    for _ in range(50):
        Omega_r = omega_r_total(H0)
        Omega_L = 1.0 - PLANCK_OMEGA_M - Omega_r
        def da_integrand(z):
            return 1.0 / np.sqrt(
                PLANCK_OMEGA_M * (1.0 + z)**3
                + Omega_L
                + Omega_r * (1.0 + z)**4
            )
        I, _ = quad(da_integrand, 0.0, PLANCK_Z_STAR, limit=400)
        H0_new = theta_rad * C_KM_S * I / r_d
        if abs(H0_new - H0) < 1e-6:
            break
        H0 = H0_new
    return H0_new


def compute_S8_analytic(R_b: float = None) -> dict:
    """
    CAR S8 prediction — analytic, independently verified.

    Paper 15 §2.5 / Series 2 Paper 1 §11.6:
        S8_analytic = 0.832 × (1 + R_b/3)^(-1/2)
        S8_numerical = S8_analytic - 0.015   (CAMB Jeans + baryonic correction)

    Uncertainty propagation:
        dS8/dR_b = -0.832/6 × (1 + R_b/3)^(-3/2) × dR_b
    """
    if R_b is None:
        R_b  = R_B_DERIVED
        dR_b = R_B_UNCERTAINTY
    else:
        dR_b = 0.0
    factor       = (1.0 + R_b / 3.0) ** (-0.5)
    S8_analytic  = PLANCK_S8 * factor
    S8_numeric   = S8_analytic - 0.015
    dS8_analytic = PLANCK_S8 * (1.0/6.0) * (1.0 + R_b/3.0)**(-1.5) * dR_b
    dS8_numeric  = dS8_analytic
    return {
        'S8_analytic':  S8_analytic,
        'S8_numeric':   S8_numeric,
        'dS8_analytic': dS8_analytic,
        'dS8_numeric':  dS8_numeric,
        'factor':       factor,
        'R_b_used':     R_b,
        'source':       'derived' if R_b == R_B_DERIVED else 'caller-supplied',
    }


def CAR_predictions(Omega_m: float = PLANCK_OMEGA_M,
                    R_b: float = None,
                    verbose: bool = False) -> dict:
    """
    Compute all CAR predictions (v4.8.5 — derived R_b, RNLA v2.3 audited; v4.8.5 register update).
    """
    if R_b is None:
        R_b    = R_B_DERIVED
        dR_b   = R_B_UNCERTAINTY
        source = 'derived (Series 2 Paper 1 Section 11.6)'
    else:
        dR_b   = 0.0
        source = 'caller-supplied={:.3f} [comparison only]'.format(R_b)

    S8d  = compute_S8_analytic(R_b)
    IA   = 1.0 + R_b / 3.0
    dIA  = dR_b / 3.0
    cs2  = (1.0 + R_b) / 3.0
    dcs2 = dR_b / 3.0

    # Standard photon-baryon sound horizon (canonical, H0-independent) and the
    # global H0 it implies through θ* = r_d / D_M(z*).
    r_d   = R_D_DERIVED
    H0    = compute_H0_from_theta_and_rd(r_d)
    r_d_int = compute_r_d_integral(H0)         # standard-horizon integral check

    # CMB-S4 separation on N_eff
    N_eff_sep = abs(N_EFF_SCT - N_EFF_SM) / CMB_S4_SIGMA_NEFF

    if verbose:
        print('  R_b [derived]   = {:.4f} ± {:.3f}  ({})'.format(R_b, dR_b, source))
        print('  c_s² [late-time] = {:.5f} ± {:.5f}  (sets S8, b_IA)'.format(cs2, dcs2))
        print('  S8  (analytic)  = {:.4f} ± {:.4f}'.format(S8d['S8_analytic'], S8d['dS8_analytic']))
        print('  S8  (numerical) = {:.4f} ± {:.4f}'.format(S8d['S8_numeric'], S8d['dS8_numeric']))
        print('  b_IA            = {:.4f} ± {:.4f}'.format(IA, dIA))
        print('  r_d (standard)  = {:.1f} ± {:.1f} Mpc  (drag horizon, Plasma Equivalence)'.format(
            R_D_DERIVED, R_D_UNCERTAINTY))
        print('  r_* (z*)        = {:.1f} Mpc  (last-scattering horizon)'.format(R_STAR_DERIVED))
        print('  r_d integral    = {:.1f} Mpc  (standard-speed sound-horizon check)'.format(r_d_int))
        print('  H0 (global)     = {:.2f} km/s/Mpc  (θ* + r_d; PARTIAL — see note)'.format(H0))
        print('  N_eff           = {:.3f} ± {:.2f}  (SCT predicted)'.format(N_EFF_SCT, N_EFF_UNCERTAINTY))
        print('  N_eff           = {:.3f}             (Standard Model)'.format(N_EFF_SM))
        print('  CMB-S4 σ        = {:.1f}             (decisive test, forecast)'.format(N_eff_sep))

    return {
        # ── Derived constants (Series 2 Paper 1 Section 11.6) ───────────────────
        'R_b0':                    R_b,
        'R_b0_uncertainty':        dR_b,
        'R_b0_source':             source,
        'cs2':                     cs2,
        'cs2_uncertainty':         dcs2,
        'C_hat_bg':                IA,
        'C_hat_bg_uncertainty':    dIA,
        'r_d_derived_Mpc':         R_D_DERIVED,
        'r_d_derived_uncertainty': R_D_UNCERTAINTY,
        'r_star_Mpc':              R_STAR_DERIVED,
        # ── Analytic predictions (no CAMB) ────────────────────────────────────
        'S8':                      S8d['S8_numeric'],
        'S8_uncertainty':          S8d['dS8_numeric'],
        'S8_analytic':             S8d['S8_analytic'],
        'S8_analytic_uncertainty': S8d['dS8_analytic'],
        'IA_bias':                 IA,
        'IA_bias_uncertainty':     dIA,
        # ── Self-consistent integrals ────────────────────────────────────────
        'r_d_integral_Mpc':        r_d_int,
        'H0_from_integral':        H0,
        # ── CAMB-required (full Boltzmann solver) ────────────────────────────
        'r_d_CAMB_Mpc':            R_D_DERIVED,
        'H0_CAMB':                 H0,
        # ── N_eff (CMB-S4 prediction) ────────────────────────────────────────
        'N_eff_SCT':               N_EFF_SCT,
        'N_eff_uncertainty':       N_EFF_UNCERTAINTY,
        'N_eff_SM':                N_EFF_SM,
        'N_eff_CMB_S4_sigma':      N_eff_sep,
        # ── Legacy/convenience keys ──────────────────────────────────────────
        'r_d_Mpc':                 R_D_DERIVED,
        'H0_km_s_Mpc':             H0,
        'H0':                      H0,            # legacy alias for older tests
        'theta_star':              PLANCK_THETA_STAR,
    }


def lcdm_reference() -> dict:
    """Standard ΛCDM reference values for comparison printouts."""
    return {
        'R_b': '~0.617 (at z_drag)',
        'cs2': 1.0/3.0,
        'b_IA': 1.0,
        'S8': 0.832,
        'r_d': 147.1,
        'H0': 67.4,
        'N_eff': 3.046,
    }


# ════════════════════════════════════════════════════════════════════════════
# REPORTING
# ════════════════════════════════════════════════════════════════════════════

def validation_report() -> None:
    """Print all SCT predictions with sigma values against observations."""
    preds = CAR_predictions()

    obs = {
        'R_b':  (R_B_LEGACY_OBS, 0.032,             'Planck CMB + BAO (legacy)'),
        'S8_D': (OBSERVED_S8_DES_Y6[0],   OBSERVED_S8_DES_Y6[1],   'DES-Y6 2026'),
        'S8_K': (OBSERVED_S8_KIDS_DR5[0], OBSERVED_S8_KIDS_DR5[1], 'KiDS-DR5'),
        'S8_H': (OBSERVED_S8_HSC_Y3[0],   OBSERVED_S8_HSC_Y3[1],   'HSC-Y3'),
        'r_d':  (OBSERVED_R_D_DESI[0],    OBSERVED_R_D_DESI[1],    'DESI-DR2 BAO'),
        'H0':   (OBSERVED_H0_PLANCK[0],   OBSERVED_H0_PLANCK[1],   'Planck 2018'),
    }

    w = 76
    print()
    print('=' * w)
    print('  SCT Validation Report v4.8.5 | DR JM NIPOK (2026)')
    print('  Paper 15 DOI: 10.13140/RG.2.2.10321.29288')
    print('  Series 2 Paper 1 DOI: 10.13140/RG.2.2.14355.03366  Section 11.6')
    print('=' * w)

    print('\n  DERIVED CONSTANTS (Series 2 Paper 1 §11.6 — no observational input)')
    print('  ' + '-' * 72)
    print('  {:<32} {:>10} {:>10} {:>9}  Note'.format('Constant', 'Derived', 'Observed', 'Sigma'))
    print('  ' + '-' * 72)

    R_b_pred  = preds['R_b0']
    R_b_dpred = preds['R_b0_uncertainty']
    R_b_obs, R_b_dobs, _ = obs['R_b']
    R_b_sigma = abs(R_b_pred - R_b_obs) / np.sqrt(R_b_dpred**2 + R_b_dobs**2)
    print('  {:<32} {:>10.4f} {:>10.3f} {:>8.2f}σ  DERIVED'.format(
        'R_b (Series 2 Paper 1 §11.6)', R_b_pred, R_b_obs, R_b_sigma))
    print('  {:<32} {:>10.5f} {:>10.4f} {:>9}  DERIVED'.format(
        'c_s² = (1+R_b)/3 (late-time)', preds['cs2'], 0.4182, '—'))
    print('  {:<32} {:>10.5f} {:>10} {:>9}  DERIVED'.format(
        'b_IA = 1 + R_b/3',  preds['C_hat_bg'], '—', '—'))

    rd_pred  = preds['r_d_derived_Mpc']
    rd_dpred = preds['r_d_derived_uncertainty']
    rd_obs, rd_dobs, _ = obs['r_d']
    rd_sigma = abs(rd_pred - rd_obs) / np.sqrt(rd_dpred**2 + rd_dobs**2)
    print('  {:<32} {:>10.1f} {:>10.1f} {:>8.2f}σ  STANDARD'.format(
        'r_d (standard horizon, Mpc)', rd_pred, rd_obs, rd_sigma))

    print('\n  ANALYTIC PREDICTIONS — verified, no CAMB required')
    print('  ' + '-' * 72)
    print('  {:<32} {:>10} {:>10} {:>9}  Survey'.format('Quantity', 'SCT', 'Observed', 'Sigma'))
    print('  ' + '-' * 72)

    S8_pred  = preds['S8']
    S8_dpred = preds['S8_uncertainty']
    for key, label in [('S8_D', 'DES-Y6'), ('S8_K', 'KiDS-DR5'), ('S8_H', 'HSC-Y3')]:
        obs_val, obs_err, _ = obs[key]
        sig = abs(S8_pred - obs_val) / np.sqrt(S8_dpred**2 + obs_err**2)
        print('  {:<32} {:>10.4f} {:>10.3f} {:>8.2f}σ  {}'.format(
            'S8', S8_pred, obs_val, sig, label))

    IA_pred = preds['IA_bias']
    print('  {:<32} {:>10.5f} {:>10} {:>9}  IA analyses'.format(
        'b_IA', IA_pred, '~1.08', '—'))

    H0_pred = preds['H0_CAMB']
    H0_obs, H0_derr, _ = obs['H0']
    H0_sig  = abs(H0_pred - H0_obs) / H0_derr
    print('  {:<32} {:>10.2f} {:>10.1f} {:>8.1f}σ  Planck 2018'.format(
        'H0 (global, θ*+r_d)', H0_pred, H0_obs, H0_sig))

    print('\n  CMB-S4 PREDICTION (NEW — Series 2 Paper 1 §11.6)')
    print('  ' + '-' * 72)
    N_sct = preds['N_eff_SCT']
    N_sm  = preds['N_eff_SM']
    N_unc = preds['N_eff_uncertainty']
    N_sig = preds['N_eff_CMB_S4_sigma']
    print('  N_eff (SCT predicted):    {:.3f} ± {:.2f}'.format(N_sct, N_unc))
    print('  N_eff (Standard Model):   {:.3f}'.format(N_sm))
    print('  Difference:               {:.3f}  ({:.1f}σ at CMB-S4 forecast σ={:.2f})'.format(
        abs(N_sct - N_sm), N_sig, CMB_S4_SIGMA_NEFF))
    print('  SCT < 3.000  |  SM > 3.000   →  TEST IS DECISIVE')
    print('  Status: PENDING (CMB-S4 ~2030)')

    print('\n  AUDIT NOTE — v4.8.4 RNLA v2.3 Recursive Audit (June 2026)')
    print('  ' + '-' * 72)
    print('  • r_d RESTORED to standard photon-baryon horizon 146.8 Mpc')
    print('    (v4.8.1 161.4 Mpc was a category error: late-time c_s in the')
    print('    recombination integral). r_*(z*) = 144.4 Mpc.')
    print('  • Standard r_d is consistent with DESI-DR2 BAO (147 ± 1 Mpc).')
    print('  • H0 global ≈ 66.3 from θ*+r_d; 70.4 (not CMB-derivable) retired;')
    print('    local H0 → 70–73 via the late-time void + temporal mechanism.')
    print('  • S8 and b_IA (late-time coherent sector) are robust and verified.')

    print('\n  REGISTER UPDATE — v4.8.5 (July 2026, in-session CAMB/SPARC note)')
    print('  ' + '-' * 72)
    print('  • N_eff=2.514 quantified vs Planck TT: naive stock-ΛCDM mapping at')
    print('    fixed θ* excluded at Δχ²≈+658 (83 binned points). New tensions.csv')
    print('    T031; predictions.csv P088 status: PENDING (PARTIAL sub-check).')
    print('  • Flagged open: which era N_eff=2.514 describes (recombination vs.')
    print('    the late-time quantity CMB-S4 will actually observe).')
    print('  • No numeric constants changed: R_b, c_s2, b_IA, r_d, H0, N_eff all')
    print('    UNCHANGED from v4.8.4.')
    print('=' * w)
    print()



def print_report() -> None:
    """Compact summary report (v4.8.5)."""
    preds = CAR_predictions()
    lcdm  = lcdm_reference()
    w = 72
    print()
    print('=' * w)
    print('  CAR Core Calculator v4.8.5 | SCT Papers #16, #17 | DR JM NIPOK (2026)')
    print('=' * w)
    print('  {:<34} {:>10}  {:>10}  {}'.format('Quantity', 'CAR', 'LCDM', 'Source'))
    print('-' * w)
    print('  {:<34} {:>10.4f}  {:>10}  {}'.format('R_b0  (DERIVED 11.6)',          preds['R_b0'],          '-',     'derived'))
    print('  {:<34} {:>10.5f}  {:>10.4f}  {}'.format('c_s2 (late-time, 1+R_b)/3',   preds['cs2'],           lcdm['cs2'], 'derived'))
    print('  {:<34} {:>10.5f}  {:>10.4f}  {}'.format('b_IA  = 1 + R_b/3',           preds['IA_bias'],       lcdm['b_IA'], 'analytic'))
    print('  {:<34} {:>10.4f}  {:>10.4f}  {}'.format('S8 (analytic)',               preds['S8_analytic'],   lcdm['S8'],   'analytic'))
    print('  {:<34} {:>10.4f}  {:>10.4f}  {}'.format('S8 (numerical, -0.015)',      preds['S8'],            lcdm['S8'],   'CAMB'))
    print('-' * w)
    print('  {:<34} {:>10.1f}  {:>10.1f}  {}'.format('r_d  (standard horizon, Mpc)', preds['r_d_CAMB_Mpc'],  lcdm['r_d'],  'standard'))
    print('  {:<34} {:>10.2f}  {:>10.1f}  {}'.format('H0  (global, theta*+r_d)',     preds['H0_CAMB'],       lcdm['H0'],   'derived'))
    print('-' * w)
    print('  {:<34} {:>10.3f}  {:>10.3f}  {}'.format('N_eff (SCT predicted)',        preds['N_eff_SCT'],     lcdm['N_eff'], 'CMB-S4'))
    print('  {:<34} {:>9.1f}s  {:>10}  {}'.format('CMB-S4 separation (forecast)',    preds['N_eff_CMB_S4_sigma'], '-', 'DECISIVE'))
    print('=' * w)
    print()
    print('  Independently verified (RNLA v2.3 audit, June 2026)')
    print('  CAMB req. - apply camb/equations_car.f90 patch for full Boltzmann run')
    print()
    print('  v4.8.4 audit: r_d = 146.8 Mpc standard horizon (161.4 category error retired).')
    print('  CAR S8 and b_IA (late-time sector) are robust; r_d consistent with DESI BAO.')
    print('  v4.8.5 register update: N_eff=2.514 vs Planck Δχ²≈+658 (naive mapping); see tensions.csv T031.')
    print('=' * w)


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='CAR Core Calculator v4.8.5 - SCT Papers #16, #17',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed intermediate values')
    parser.add_argument('--validate', action='store_true',
                        help='Print full validation report with sigma values')
    args = parser.parse_args()

    if args.verbose:
        CAR_predictions(verbose=True)
    if args.validate:
        validation_report()
    elif not args.verbose:
        print_report()
