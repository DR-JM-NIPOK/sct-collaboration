# SCT Cosmology Series
## Codified Acoustic Relation (CAR) — Parameter-Free Cosmology

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.13140%2FRG.2.2.10321.29288-blue)](https://doi.org/10.13140/RG.2.2.10321.29288)
[![OSF](https://img.shields.io/badge/OSF-t8zny-green)](https://osf.io/t8zny)

This repository contains all code, configuration, and data-download scripts for the
**Successive Collision Theory (SCT) Cosmology Series**, culminating in the
**Codified Acoustic Relation (CAR)** framework presented in Paper #16.

---

### Paper Reference

NIPOK, DR JM, *"From Chaos to Codified Acoustics: A Parameter-Free Collision Geometry
That Unifies DESI-DR2, DES-Y6, HSC-Y3, and KiDS-DR5 While Resolving H₀ and S₈ Tensions,"*
SCT Cosmology Series Paper #16, N.J.I.T. (2026).

- **DOI**: [10.13140/RG.2.2.10321.29288](https://doi.org/10.13140/RG.2.2.10321.29288)
- **OSF Archive**: [https://osf.io/t8zny](https://osf.io/t8zny)
- **Preprint**: [thenaturalstateofnature.org/PREPRINTS/](https://thenaturalstateofnature.org/PREPRINTS/)

---

### Theory in One Equation

The CAR master equation for the primordial plasma sound speed:

```
c_s²(z) = [1 + R_b(z)] / 3       R_b(z) = R_b0 / (1+z)
```

where `R_b0 = 0.2545 ± 0.032` is DERIVED from SO(3) cascade geometry and QCD
Israel-Darmois junction conditions (Series 2 Paper 1 Section 11.6 — no observational input).
DO NOT compute R_b0 from Ω_b h² / Ω_γ h².

This single relation — **zero free parameters** — predicts:

| Quantity | CAR Prediction | Observation | Tension |
|---|---|---|---|
| r_d (standard) | 146.8 ± 5 Mpc | DESI-DR2: 147.0 ± 1.0 Mpc | 0.2σ |
| H₀ (global) | 66.3 km/s/Mpc | Planck: 67.4 ± 0.5 | ~2σ; local→70-73 via void+temporal |
| S₈ | 0.783 ± 0.015 | DES-Y6: 0.780 ± 0.012 | <0.2σ |
| b_IA | 1.0848 ± 0.011 | DES-Y6 fit: 1.08 ± 0.04 | <0.2σ |

---

### Repository Structure

```
sct-collaboration/
├── sct_core.py              ← Core CAR calculator (run this first)
├── requirements.txt         ← Python dependencies
├── README.md                ← This file
├── SETUP_INSTRUCTIONS.md    ← Step-by-step environment setup
├── LICENSE                  ← GPL-3.0
│
├── camb/                    ← Modified CAMB Boltzmann solver
│   ├── README.md
│   ├── equations_CAR.patch  ← Unified diff for equations.f90
│   └── equations_car_test.py
│
├── class/                   ← Modified CLASS Boltzmann solver
│   ├── README.md
│   ├── perturbations_CAR.patch  ← Unified diff for perturbations.c
│   └── class_car_test.py
│
├── likelihoods/             ← Likelihood modules for each survey
│   ├── README.md
│   ├── car_likelihood_base.py   ← Abstract base class
│   ├── desi_dr2_bao.py          ← DESI-DR2 BAO likelihood
│   ├── des_y6_3x2pt.py          ← DES-Y6 3×2pt likelihood
│   ├── hsc_y3_wl.py             ← HSC-Y3 weak lensing likelihood
│   ├── kids_dr5_wl.py           ← KiDS-DR5 weak lensing likelihood
│   └── planck_pr4_cmb.py        ← Planck PR4 CMB compressed likelihood
│
├── chains/                  ← PolyChord nested sampling
│   ├── README.md
│   ├── run_chains.sh            ← Master run script
│   ├── config_car.ini           ← CAR PolyChord settings
│   ├── config_lcdm.ini          ← ΛCDM PolyChord settings (comparison)
│   ├── compute_evidence.py      ← Log-evidence extraction and comparison
│   └── plot_posteriors.py       ← GetDist posterior plots
│
├── figures/                 ← All paper figure scripts
│   ├── README.md
│   ├── make_all_figures.py      ← Master figure script
│   ├── fig1_tension_summary.py  ← H₀, S₈, r_d tension plot
│   ├── fig2_sound_horizon.py    ← c_s(z) and r_d integral
│   ├── fig3_posteriors.py       ← 2D posterior contours
│   ├── fig4_residuals.py        ← Standardised residuals by dataset
│   └── fig5_forecasts.py        ← Euclid/LSST/CMB-S4 forecasts
│
├── data/                    ← Data vectors and download scripts
│   ├── README.md
│   ├── download_all.sh          ← Fetch all public data
│   ├── desi_dr2_bao_data.py     ← DESI-DR2 BAO data loader
│   ├── des_y6_data.py           ← DES-Y6 data loader
│   └── mock_data_generator.py   ← Generate mock data for testing
│
├── docker/                  ← Container for full reproducibility
│   ├── Dockerfile.txt
│   └── docker-compose.yml
│
└── tests/                   ← Unit and integration tests
    ├── test_sct_core.py
    ├── test_likelihoods.py
    └── test_sound_horizon.py
```

---

### Quick Start

#### Option 1: Core Calculator Only (no external dependencies)

```bash
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration
pip install numpy scipy
python sct_core.py
```

Expected output:
```
========================================================
  Codified Acoustic Relation (CAR) — Predictions
========================================================
  Quantity               CAR          ΛCDM         Δ
--------------------------------------------------------
  R_b0                  0.2545         —      (derived, Series 2 Paper 1)
  c_s²(z=0)  [×c²]      0.41817      0.33333  +0.08484  (late-time; R_b=0.2545)
  r_d  [Mpc]           146.80       147.10     -0.30   (standard horizon)
  H₀  [km/s/Mpc]        66.30        67.40     -1.10   (global; θ*+r_d)
  S₈  (numerical)        0.783        0.832    -0.049
  b_IA                   1.0848       1.000    +0.0848
```

#### Option 2: Full Analysis (requires CAMB/CLASS/PolyChord)

```bash
# Follow SETUP_INSTRUCTIONS.md first, then:
bash chains/run_chains.sh --data combined --model car
python figures/make_all_figures.py
```

#### Option 3: Docker (complete environment)

```bash
docker build -t sct-cosmology docker/
docker run -v $(pwd)/output:/output sct-cosmology bash chains/run_chains.sh
```

---

### Installation

```bash
pip install -r requirements.txt
```

For CAMB and CLASS modifications, see `camb/README.md` and `class/README.md`.

---

### Verification

The expected χ²/dof values for the CAR model against each dataset:

| Dataset | χ²/dof | p-value |
|---|---|---|
| DESI-DR2 BAO | 0.931 | 0.57 |
| Planck PR4 CMB | 1.028 | 0.41 |
| DES-Y6 Weak Lensing | 0.982 | 0.52 |
| HSC-Y3 + KiDS-DR5 | 0.994 | 0.48 |
| **Combined** | **0.990** | **0.62** |

---

### License

This repository is licensed under the **GNU General Public License v3.0** (GPL-3.0).
See [LICENSE.txt](LICENSE.txt) for full terms.

The accompanying paper is licensed under **CC BY-NC-SA 4.0**.

---

### Contact

- **Author**: DR JM NIPOK
- **ORCID**: [0009-0006-3940-4450](https://orcid.org/0009-0006-3940-4450)
- **Issues**: [GitHub Issues](https://github.com/DR-JM-NIPOK/sct-collaboration/issues)
- **OSF**: [https://osf.io/t8zny](https://osf.io/t8zny)
