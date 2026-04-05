# CHANGELOG — SCT Cosmology Repository

All notable changes to the `sct-collaboration` repository are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [3.0.0] — 2026-04-04

### Paper 17 v4.0 Epistemic Upgrade

**Reference:** Paper 17 v4.0 — DOI: [10.13140/RG.2.2.14355.03366](https://doi.org/10.13140/RG.2.2.14355.03366) — Section 11.6

**Scientific rationale:**
R_b transitions from a **matched observational parameter** (v2.0, Paper 16) to a
**derived constant of the SCT framework** (v3.0, Paper 17 v4.0 §11.6). The
derivation uses two purely geometric/field-theoretic inputs:

1. **SO(3) angular momentum structure** of the collision cascade → N_cascade = 3
   (three independent cascade planes from the SO(3) symmetry group)
2. **QCD phase transition boundary correction** from Paper 14 Israel-Darmois
   junction conditions → 13.6% energy loss at the QCD boundary

No observational input is required. The result R_b = 0.257 ± 0.032 agrees
with the previously matched value R_b = 0.260 at **0.11 sigma**, closing the
circularity in the CAR Bayes factor calculation.

**New prediction:** N_eff_SCT = 2.57 ± 0.05 (vs Standard Model 3.046). These
values are on opposite sides of 3.000. CMB-S4 will distinguish them at **16 sigma**
— the test is decisive and falsifiable.

---

### Files Changed

### Branch-specific changes — `chains`

Updated: sct_core.py, SETUP_INSTRUCTIONS.md, README.md, config_car.ini (R_b removed from prior), run_polychord.py (no Omega_b_h2 reconstruction), plot_posteriors.py

All other files on this branch are preserved exactly as originally committed.

---

#### `sct_core.py`

- **Epistemic upgrade v2.0 → v3.0** documented in module docstring
- Added derived constants block at top of file:
  - `R_B_DERIVED = 0.257` (was hardcoded `R_B0 = 0.260`)
  - `R_B_UNCERTAINTY = 0.032`
  - `CS2_DERIVED = 0.41900` (updated from 0.420; follows from derived R_b)
  - `CS2_UNCERTAINTY = 0.011`
  - `C_HAT_BG_DERIVED = 1.08567` (updated from 1.087)
  - `C_HAT_BG_UNCERTAINTY = 0.011`
  - `R_D_DERIVED = 146.8` Mpc (Paper 17 v4.0)
  - `R_D_UNCERTAINTY = 5.0` Mpc
  - `N_EFF_SCT = 2.566` (new prediction)
  - `N_EFF_UNCERTAINTY = 0.05`
  - `N_EFF_SM = 3.046`
  - `R_B_LEGACY_OBS = 0.260` (reference only — DO NOT USE AS INPUT)
- `R_B0` alias set to `R_B_DERIVED` — all internal computation uses derived value
- `CAR_predictions()`: now defaults to derived R_b; accepts optional `R_b` argument
  for comparison purposes only; never computes R_b from Omega_b_h2
- `compute_S8_analytic()`: updated to use derived R_b with full uncertainty propagation
- Added `validation_report()` function: prints all predictions with sigma values vs
  observations, N_eff CMB-S4 status, and legacy comparison table
- `print_report()`: updated to v3.0, shows derived constants and N_eff prediction
- All docstrings updated with Paper 17 v4.0 §11.6 citation
- Warning added: `R_b = 0.260 must NOT be used as input`

**Values changed:**
  - `R_B0`: 0.260 → **0.257** (derived)
  - `c_s²/c²`: 0.420 → **0.419** (derived; follows from R_b)
  - `C_hat_bg`: 1.087 → **1.086** (derived; follows from R_b)
  - `r_d` (CAMB): 149.2 → **146.8 ± 5 Mpc** (Paper 17 v4.0)
  - `N_eff_SCT`: — → **2.57 ± 0.05** (new prediction)

#### `README.md`

- Added **⚠ WARNING** block: R_b = 0.260 must not be used as input
- Added **Paper 17 v4.0** to Paper Reference section with full DOI
- Added **Key Constants table** (all DERIVED):
  - R_b = 0.257 ± 0.032 DERIVED
  - c_s² = 0.419 ± 0.011 DERIVED
  - Ĉ_bg = 1.086 ± 0.011 DERIVED
  - r_d = 146.8 ± 5 Mpc DERIVED
  - N_eff SCT = 2.57 ± 0.05 PREDICTED — CMB-S4 16σ
- Updated Verified Outputs table: R_b0 row changed from `0.260 / Matched input`
  to `0.257 ± 0.032 / DERIVED`
- Added `python sct_core.py --validate` to Quick Start
- Added `CHANGELOG.md` to Repository Contents table
- Version header updated to v3.0

#### `SETUP_INSTRUCTIONS.md`

- Added **⚠ WARNING** block at top: R_b = 0.260 must not be used as input
- Version header updated from v2.0 to v3.0
- Expected output block: `R_b0` row updated to `0.2570 (derived)`, label changed
  from `analytic` to `derived`
- Bug fix note [1]: updated from `R_b0 = 0.260` to `R_b0 = 0.257 derived
  (Paper 17 v4.0 §11.6) — DO NOT compute from Omega_b_h2`
- **Key Parameters table** (§ 9): R_b0 row updated:
  - Value: `0.260` → `0.257 ± 0.032`
  - Status: `Matched observational input` → `Derived constant — Paper 17 v4.0 §11.6`
- Troubleshooting: updated R_b0 bug description to reference derived value
- All remaining instances of `0.420` replaced with `0.419`
- All remaining instances of `0.260` replaced with `0.257 (derived, Paper 17 v4.0)`

#### `growth.py`

- Import updated: `R_B0` now sourced from `sct_core.R_B_DERIVED` (0.257),
  not from `coherence.R_B0` (which was 0.260)
- Comment added: `DO NOT use 0.260 — use R_B_DERIVED = 0.257`
- `cs2_CAR_background()` docstring: `cs² = 0.420` → `cs² = 0.419 (derived, R_b=0.257)`
- `cs2_CAR_perturbation()` docstring: `cs² = 0.420` → `cs² = 0.419 (derived, Paper 17 v4.0 §11.6)`
- Prose: `Full R_b0 = 0.260 applies` → `Full R_b0 = 0.257 (derived, Paper 17 v4.0 §11.6) applies`

**Values changed:**
  - All `cs² = 0.420` references → **0.419**
  - `R_B0` source: coherence.py (0.260) → sct_core.py derived (0.257)

#### `coherence.py`

- Line 54: `R_B0 = 0.260` → `R_B0 = 0.257` with updated comment:
  `# CAR coherence parameter (DERIVED, Paper 17 v4.0 Section 11.6)`
- Warning comment added: `DO NOT use 0.260 — that was the legacy matched value`

**Values changed:**
  - `R_B0`: 0.260 → **0.257** (derived)
  - All downstream quantities (A_eff, C_hat_background, etc.) update automatically

#### `CHANGELOG.md` (this file)

- Created new file (did not previously exist)

---

## [2.0.0] — 2026-04-01

### Three Critical Bug Fixes

**Reference:** Paper 16 — DOI: 10.13140/RG.2.2.10321.29288

#### Bug 1 — R_b0 convention (CRITICAL)
- **Broken:** `R_b0 = 4×Ω_b_h²/(3×Ω_γ_h²) = 1196.9` (physical z=0 ratio, wrong formula)
- **Fixed:** `R_b0 = 0.260` (CAR coherence parameter, Paper 16 §2.1)
- **Effect:** S8: 0.042 → 0.783; b_IA: 400 → 1.087

#### Bug 2 — θ* unit conversion (CRITICAL)
- **Broken:** `theta_star_rad = theta_star × π/180` (treats as degrees)
- **Fixed:** `theta_star_rad = theta_star / 100` (Planck reports 100×θ*)
- **Effect:** H0: 1.3 → 70.4 km/s/Mpc

#### Bug 3 — r_d integral normalisation
- **Broken:** `r_d = integral × c/100` (assumes H0=100)
- **Fixed:** `r_d = integral × c/H0`
- **Effect:** r_d: 133 → 149 Mpc

---

## [1.0.0] — 2026-03-01

### Initial release

- `sct_core.py` v1.0 — CAR Core Calculator
- Paper 16 submitted: DOI 10.13140/RG.2.2.10321.29288
