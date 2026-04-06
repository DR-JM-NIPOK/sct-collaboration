# SCT Cosmology Series
## Codified Acoustic Relation (CAR) — Derived-Parameter Cosmology

**Version 3.0** | April 2026 | Paper 17 v4.8 epistemic upgrade: R_b derived from first principles
If you downloaded this code before April 2026, please re-clone.

---

### ⚠ WARNING — DO NOT USE R_b = 0.260 AS INPUT

R_b = 0.260 was the **matched observational value** in v1.0/v2.0 (Paper 16).
As of v3.0, R_b = **0.2545 ± 0.032 is a DERIVED constant** (Paper 17 v4.8 §11.6).
Passing 0.260 to `sct_core.py` introduces a circularity that the Paper 17
derivation eliminates. The module will accept it as a legacy comparison argument
only — see `CAR_predictions(R_b=R_B_LEGACY_OBS)` for that use case.

---

### Paper References

**Paper 16:**
NIPOK, DR JM — *"From Chaos to Codified Acoustics: A Parameter-Free Collision
Geometry That Unifies DESI DR2, DES Y6, HSC Y3, and KiDS DR5 While Resolving
H₀ and S₈ Tensions"* — SCT Cosmology Series Paper #16 (2026)
DOI: [10.13140/RG.2.2.10321.29288](https://doi.org/10.13140/RG.2.2.10321.29288)

**Paper 17 v4.8 (current):**
NIPOK, DR JM — *"R_b Derived from First Principles: SO(3) Cascade Geometry and
QCD Junction Conditions Eliminate the Last Matched Parameter of CAR"* —
SCT Cosmology Series Paper #17 v4.0 (2026)
DOI: [10.13140/RG.2.2.14355.03366](https://doi.org/10.13140/RG.2.2.14355.03366)
Section 11.6 contains the derivation of R_b = 0.2545 ± 0.032.

OSF Archive: [doi.org/10.17605/OSF.IO/T8ZNY](https://doi.org/10.17605/OSF.IO/T8ZNY)
ORCID: [0009-0006-3940-4450](https://orcid.org/0009-0006-3940-4450)

---

### Quick Start

```bash
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration
pip install -r requirements.txt
python sct_core.py
python sct_core.py --validate   # full sigma comparison vs observations
```

---

### Key Constants (v3.0, all DERIVED — Paper 17 v4.8 §11.6)

| Constant | Value | Status | Source |
|----------|-------|--------|--------|
| R_b | **0.2545 ± 0.032** | DERIVED | SO(3) cascade + QCD junction conditions |
| c_s² / c² | **0.4182 ± 0.011** | DERIVED | (1 + R_b) / 3 with derived R_b |
| Ĉ_bg | **1.0848 ± 0.011** | DERIVED | 1 + R_b/3 with derived R_b |
| r_d | **146.8 ± 5 Mpc** | DERIVED | CAMB + CAR patch, Paper 17 v4.8 |
| N_eff (SCT) | **2.514 ± 0.05** | PREDICTED — CMB-S4 forecast at 17.7σ | Paper 17 v4.8 §11.6 |
| N_eff (SM) | 3.046 | Standard Model | Mangano et al. 2005 |

N_eff_SCT = 2.514 and N_eff_SM = 3.046 are on **opposite sides of 3.000**.
CMB-S4 will separate these at **17.7 sigma (forecast)** — the test is decisive.

---

### Verified Outputs (v3.0)

| Quantity | CAR | ΛCDM | Source |
|----------|-----|------|--------|
| R_b0 | **0.2545 ± 0.032** | — | DERIVED (Paper 17 v4.8 §11.6) |
| c_s² | **0.4182 ± 0.011** | 0.333 | DERIVED |
| S₈ | **0.783** | 0.832 | Analytic ✓ verified |
| b_IA | **1.0848** | 1.000 | Analytic ✓ verified |
| r_d | **146.8 ± 5 Mpc** | 147.1 Mpc | CAMB required |
| H₀ | **70.4 km/s/Mpc** | 67.4 | CAMB required |
| N_eff | **2.514 ± 0.05** | 3.046 | PREDICTED — CMB-S4 |

S₈ and b_IA are independently verified analytically.
r_d and H₀ require CAMB with `camb/equations_car.f90` patch applied.

---

### What Changed

#### v3.0 — Paper 17 v4.8 Epistemic Upgrade (April 2026)

R_b transitions from a **matched observational parameter** to a **derived constant**:

| Item | v2.0 (matched) | v3.0 (derived) | Scientific reason |
|------|----------------|----------------|-------------------|
| R_b | 0.260 (observed) | **0.2545 ± 0.032** | SO(3) cascade + QCD junction |
| c_s² | 0.420 (ΛCDM ref) | **0.4182 ± 0.011** | Follows from derived R_b |
| Ĉ_bg | 1.087 (ΛCDM ref) | **1.0848 ± 0.011** | Follows from derived R_b |
| N_eff | — | **2.514 ± 0.05 (new)** | CMB-S4 test at 16σ |

The derivation uses N_cascade = 3 from SO(3) angular momentum structure and
a 13.6% energy loss from Paper 14 Israel-Darmois junction conditions. No
observational input is required. Agreement with observed R_b = 0.260 is at
0.11σ — a post-diction that closes the Bayesian circularity.

#### v2.0 — Three Critical Bugs Corrected (April 2026)

| Bug | Broken | Fixed | Effect |
|-----|--------|-------|--------|
| R_b0 convention | 4Ω_bh²/3Ω_γh² = **1196.9** | **0.2545** (derived, P17 v4.8) | S8: 0.042→0.783 |
| θ* units | 1.04105 × π/180 (degrees) | 1.04105 / **100** (radians) | H0: 1.3→70.4 |
| r_d normalisation | integral × c/**100** | integral × c/**H0** | r_d: 133→149 Mpc |

---

### Repository Contents

| File/Directory | Description |
|---|---|
| `sct_core.py` | Core CAR calculator — **v3.0, R_b derived, all bugs fixed** |
| `camb/equations_car.f90` | CAMB Fortran patch for CAR sound speed |
| `class/perturbations_car.c` | CLASS C patch for CAR sound speed |
| `likelihoods/` | Likelihood modules for DESI, DES, HSC, KiDS, Planck |
| `chains/` | PolyChord nested sampling configuration |
| `figures/` | Scripts to reproduce all paper figures |
| `data/download_data.sh` | Script to download real survey data |
| `docker/Dockerfile` | Container for one-command reproducibility |
| `CHANGELOG.md` | Full history of all parameter and code changes |

---

### Important Notes for Independent Verification

**Real data required for valid Bayesian evidence:**
Without files in `data/`, likelihood modules use mock data. The 44:1 Bayesian
odds figure in Paper 16 requires real DESI-DR2, DES-Y6, HSC-Y3, and KiDS-DR5
data vectors — download with `data/download_data.sh`.

**r_d from simple integral vs CAMB:**
The Python integral in `sct_core.py` gives r_d ≈ 158 Mpc. The Paper 17 v4.8
derived value of 146.8 ± 5 Mpc comes from the full CAMB Boltzmann solver with
the CAR patch applied. The difference comes from tight-coupling and
diffusion-damping corrections in CAMB's full hierarchy.

---

### License

GPL-3.0 | Contact: DR JM NIPOK via GitHub Issues
