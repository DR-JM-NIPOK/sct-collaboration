# CHANGELOG — SCT Cosmology Repository

All notable changes to the `sct-collaboration` repository are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [4.8.1] — 2026-04-26

### NLA Recursive Audit

This release is the result of a comprehensive Numerical, Logical, and
Attributional (NLA) recursive audit of every file in the repository.
"Recursive" here means: when one canonical value changed, every dependent
claim was traced and updated. The audit identified **eleven inconsistencies**
across files that previously claimed to be at "v4.8". They are all resolved
in this release.

### Audit Findings — what was wrong before v4.8.1

| # | Finding | Severity | Status |
|---|---------|----------|--------|
| 1 | Two different `sct_core.py` files in repo (root v4.8 derived vs. branches still on v2.0 matched-R_b) | CRITICAL | FIXED |
| 2 | `camb/equations_car.f90` Fortran patch used standard density-ratio R (≈ 673 at z=0) inside the CAR formula, producing superluminal cs² and unphysical r_drag ≈ 182 Mpc | CRITICAL | FIXED |
| 3 | `class/perturbations_CAR.patch` had the same bug as #2 | CRITICAL | FIXED |
| 4 | Hardcoded `R_D_DERIVED = 146.8 Mpc` in sct_core.py was based on a CAMB run that could not be reproduced; canonical value is 161.4 Mpc | CRITICAL | FIXED |
| 5 | Hardcoded `'H0_CAMB': 70.4` returned without verification — the self-consistent θ* iteration is mathematically degenerate at fixed Ωm and converges spuriously | CRITICAL | DOCUMENTED |
| 6 | `predictions.csv` row P03 stated `r_d = 149.2 ± 0.4 Mpc` and "CONFIRMED" — does not match any reproducible CAMB or Python output | CRITICAL | FIXED |
| 7 | `predictions.csv` row P04 stated `b_IA = 1.087 ± 0.005` — should be 1.0848 ± 0.011 (derived from R_b = 0.2545) | FLAG | FIXED |
| 8 | `validation_report()` in sct_core.py used `S8_K = 0.815` for KiDS-DR5 — that was an internal preliminary value; published is 0.788 ± 0.014 | CRITICAL | FIXED |
| 9 | `equations_car_test.py` asserted `144.0 < r_d < 161.0`, a tolerance loose enough to pass for unmodified ΛCDM | FLAG | FIXED |
| 10 | `equations_car_test.py` declared "expected r_drag ≈ 149.1 Mpc" with the comment "if you get ~150 Mpc the patch was not applied" — both numbers wrong for stated parameters | CRITICAL | FIXED |
| 11 | CHANGELOG.md still described v4.0 (R_b = 0.257, N_eff = 2.57) while code was at v4.8 (R_b = 0.2545, N_eff = 2.514) | FLAG | FIXED |

### Major scientific findings of the audit

The audit established that **the CAR ansatz, when implemented honestly,
predicts r_d ≈ 161.4 Mpc**, not the previously claimed 146.8 Mpc or 149.1 Mpc.
This means the CAR ansatz **does NOT close the DESI-DR2 BAO tension**
(observed r_d = 147 ± 1 Mpc); under canonical CAR, r_d sits 14σ above DESI
and 28σ above Planck. By contrast:

- **S8 = 0.7838 ± 0.015** remains robust and matches DES-Y6 (0.780 ± 0.012),
  KiDS-DR5 (0.788 ± 0.014), and HSC-Y3 (0.776 ± 0.020) within 1σ.
- **b_IA = 1.0848 ± 0.011** matches IA analyses to better than 0.5σ.
- **N_eff_SCT = 2.514 ± 0.05** (predicted) remains a clean falsifiable
  CMB-S4 test, separated from the SM N_eff = 3.046 by 17.7σ at the
  forecast CMB-S4 sensitivity (σ_N_eff ≈ 0.03).

The CAR ansatz therefore needs reframing as resolving S8 and contributing
to (but not closing) the H0 tension, while the BAO sector requires
additional physics. **Paper 16 v3.0 will rewrite the BAO discussion**;
S8/b_IA results carry forward unchanged.

### Changed files in v4.8.1

#### `sct_core.py` (root) — synchronized to canonical
- `R_D_DERIVED`: 146.8 → **161.4** Mpc (with R_D_UNCERTAINTY 5.0 → 0.3)
- `H0_CAMB` hardcoding: removed (now derived from self-consistent integral
  iteration with documented degeneracy caveat)
- `validation_report()`: KiDS-DR5 S8 from 0.815 → **0.788 ± 0.014** (published value)
- `OBSERVED_*` constants added with provenance for transparent comparison
- Documentation rewritten with explicit r_d audit note
- Version stamps and CLI banner updated to v4.8.1

#### Branch `sct_core.py` copies — replaced
- `camb/sct_core.py`, `chains/sct_core.py`, `class/sct_core.py`,
  `data/sct_core.py`, `docker/sct_core.py`, `figures/sct_core.py`,
  `likelihoods/sct_core.py`, `tests/sct_core.py`,
  `paper17-v4-deriver-Rb/sct_core.py`
- All replaced with byte-identical copy of canonical root sct_core.py.
  Previously these were on a v2.0/v3.0 hybrid version with R_b = 0.260
  and the broken 4Ω_b/(3Ω_γ) formula.

#### `camb/equations_car.f90` — corrected
- Replaced `R = 3*grho_b/(4*grho_g); cs2 = (1+R)/3` (which gave
  superluminal cs² at z=0, unphysical) with the canonical CAR
  prescription using SCT-derived R_b(z) = R_B_DERIVED/(1+z).
- Now produces r_drag ≈ 161.4 Mpc when applied to CAMB at H0=70.4.

#### `camb/equations_CAR.patch` — corrected
- Updated patch text to reflect Fortran fix above.
- Patch now contains explicit warning: "do NOT reuse standard R in CAR formula".

#### `camb/equations_car_test.py` — corrected
- Removed loose `144.0 < r_d < 161.0` tolerance (which passed for unmodified ΛCDM).
- New tolerance: `|r_d - 161.4| < 0.5` for canonical Python check.
- Replaced erroneous "if you get ~150 Mpc the patch was not applied" comment.
- Added `test_canonical_constants()` for first-line sanity check.
- Added `test_camb_patched()` clearly distinguishing PASS/INFO/FAIL with
  diagnostic guidance for each failure mode.

#### `class/perturbations_CAR.patch` — corrected (same fix as Fortran)
#### `class/class_car_test.py` — corrected (same fix as CAMB test)

#### `predictions.csv` — corrected
- P01: S8 0.783 → **0.7838** (full precision from canonical analytic chain)
- P03: r_d 149.2 ± 0.4 "CONFIRMED" → **161.4 ± 0.3 "DOES NOT CLOSE TENSION"**
- P04: b_IA 1.087 ± 0.005 → **1.0848 ± 0.011** (canonical, audit-corrected)
- P31 (NEW): N_eff_SCT = 2.514 ± 0.05 from Paper 17 v4.8 §11.6
- P32 (NEW): R_b derived = 0.2545 ± 0.032 from Paper 17 v4.8 §11.6

#### `CHANGELOG.md` (this file) — fully synchronized to v4.8.1
#### `README.md` (root) — updated with v4.8.1 status and audit summary

### Files not changed in v4.8.1

These were checked by the audit and confirmed to be already consistent:
- `bh_interior.py`, `coherence.py`, `growth.py`, `hereditary.py` — no
  CAR-dependent values
- `tensions.csv` — observational data, independent of theoretical chain
- `data/mock_data_generator.py`, `figures/make_all_figures.py` —
  use sct_core.py, automatically inherit corrections via update
- `chains/run_polychord.py`, `chains/compute_evidence.py`,
  `chains/plot_posteriors.py` — implementation-specific, not affected by
  R_b value (they sample a likelihood already constructed from sct_core.py)
- `chains/config_*.ini` — sampling parameters, no CAR values
- `likelihoods/*.py` — call sct_core.py, automatically inherit corrections

### How to verify the audit

After cloning v4.8.1:

```bash
python3 sct_core.py --validate     # full report
python3 -m pytest tests/           # all tests should pass
python3 camb/equations_car_test.py # standalone (no CAMB) verification
```

The `--validate` output should show:
- R_b = 0.2545 ± 0.032 derived
- c_s² = 0.41817 derived
- b_IA = 1.0848 derived
- S8 = 0.7838 (numerical), 0.7988 (analytic)
- r_d = 161.40 Mpc (canonical CAR)
- N_eff = 2.514 ± 0.05 (CMB-S4 separation 17.7σ)

---

## [4.8.0] — 2026-04-06

### Epistemic upgrade — R_b transitions to derived constant

Paper 17 v4.8 Section 11.6 derives R_b = 0.2545 ± 0.032 from first
principles: SO(3) cascade angular momentum + QCD junction conditions
(Israel-Darmois, Paper 14, 13.6% loss). No observational input.

- `R_B_DERIVED = 0.2545` replaces matched value 0.260
- `c_s² = 0.41817` derived from (1 + R_b)/3
- `C_hat_bg = 1.08483` derived from 1 + R_b/3
- `N_eff_SCT = 2.514 ± 0.05` predicted (CMB-S4 testable)

⚠ Known issue at v4.8.0 (resolved in v4.8.1): the CAMB Fortran patch
and most branch copies of sct_core.py were not synchronized. v4.8.1
performs the full sweep.

---

## [4.0.0] — 2026-04-01 (superseded by v4.8.0/v4.8.1)

Initial epistemic upgrade — R_b = 0.257 derivation (later refined to 0.2545
in v4.8). N_eff prediction = 2.57 (later refined to 2.514).

---

## [2.0.0] — 2026-04 (superseded)

Three critical bugs corrected:
- Bug 1: R_b0 convention (4·Ω_b·h² / 3·Ω_γ·h² is NOT the SCT R_b)
- Bug 2: theta_star unit conversion (Planck reports 100·θ*, not degrees)
- Bug 3: r_d integral normalization (× c/H0, not × c/100)

---

## [1.0.0] — 2026-03

Original release. R_b = 0.260 as matched observational parameter.
Subsequently identified as having internal circularity in Bayesian
evidence calculations (now resolved by v4.8 derivation).

---

## Authoring & License

DR JM NIPOK | N.J.I.T. | ORCID 0009-0006-3940-4450
Paper 16 DOI: 10.13140/RG.2.2.10321.29288
Paper 17 DOI: 10.13140/RG.2.2.14355.03366
License: GPL-3.0
