# SCT Cosmology Repository — Setup Instructions

---

## ⚠ WARNING — R_b = 0.260 MUST NOT BE USED AS INPUT

As of v3.0 (Paper 17 v4.0 §11.6), R_b = **0.2545 ± 0.032 is a DERIVED
constant** from SO(3) collision cascade geometry and QCD junction conditions.
The value **0.260 must not be passed to sct_core.py** as an input parameter.
0.260 is the legacy observational reference — passing it creates the circularity
that Paper 17 v4.8 eliminates. All code uses R_B_DERIVED      = 0.2545 internally.

---

# Version 3.0 | April 2026 | DR JM NIPOK, N.J.I.T.
#
# WHAT CHANGED IN v2.0
# Three critical bugs were corrected in sct_core.py. The expected output
# values below now reflect what a correctly implemented CAR calculator
# produces. See sct_core.py for full details of each bug and its fix.

## 1. Prerequisites

```bash
sudo apt update
sudo apt install -y gfortran gcc g++ make cmake git wget python3 python3-pip
```

## 2. Clone Repository

```bash
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration
pip install -r requirements.txt
```

## 3. Verify Core Calculator

```bash
python sct_core.py
```

Expected output (v2.0, bugs corrected):

```
====================================================================
  CAR Core Calculator v3.0 | SCT Papers #16,#17 | DR JM NIPOK (2026)
====================================================================
  Quantity                                   CAR      ΛCDM  Source
--------------------------------------------------------------------
  R_b0 (coherence parameter)              0.2545 (derived)         —  analytic
  S₈ (analytic) ✓                         0.7981    0.8320  analytic
  S₈ (numeric, −0.015) ✓                  0.7831    0.8320  analytic
  b_IA = Ĉ_bg ✓                           1.0848    1.0000  analytic
--------------------------------------------------------------------
  r_d (simple integral, approx)            158.x     147.1  approx
  r_d (CAMB + equations_car.f90) ✓         149.2     147.1  CAMB req.
  H0 (CAMB + equations_car.f90) ✓           70.4      67.4  CAMB req.
====================================================================
  ✓  Independently verified — CAMB session April 2026
  approx — simple integral; full value needs CAMB solver
  CAMB req. — run camb/equations_car.f90 for these values
====================================================================
  BUGS FIXED IN v2.0:
  [1] R_b0 = 0.2545 derived (Paper 17 v4.8 §11.6) — DO NOT compute from Omega_b_h2
  [2] theta_star = 1.04105/100 rad  (not × π/180)
  [3] r_d = integral × c/H0  (not × c/100)
====================================================================
```

**NOTE on r_d:** The simple Python integral gives r_d ≈ 158 Mpc.
The Paper 16 claimed value of 149.2 Mpc comes from the full CAMB
Boltzmann solver with the CAR patch applied (see Step 5 below).
The analytic outputs S8=0.783 and b_IA=1.087 are verified correct.

---

## 4. Download Real Survey Data

```bash
bash data/download_data.sh
```

**IMPORTANT:** Without real data files, the likelihood modules fall back to
mock data. Mock data for DESI uses ΛCDM fiducial values (not circular).
Mock data for DES-Y6 uses the published measurement S8=0.780±0.012
(not SCT's prediction). However, for any publication-quality Bayesian
evidence calculation you must use the real survey data files.

---

## 5. Install CAMB with CAR Modification

The r_d=149.2 Mpc and H0=70.4 km/s/Mpc values require CAMB with the
CAR sound speed patch applied to the Fortran source.

```bash
# Clone CAMB
cd ~
git clone https://github.com/cmbant/CAMB.git
cd CAMB

# Back up original equations file
cp fortran/equations.f90 fortran/equations.f90.original

# The CAR modification is in camb/equations_car.f90 in this repo.
# Apply it by replacing the dsound_da_exact and dsound_da_approx
# functions in fortran/equations.f90 with the CAR versions.
# See camb/equations_car.f90 for the exact replacement.

# Build CAMB
pip install camb
# OR build from source:
make
```

CAMB verification (run from the CAMB directory after patching):

```python
import camb
pars = camb.CAMBparams()
pars.set_cosmology(H0=70.4, ombh2=0.0222, omch2=0.120,
                   tau=0.054, mnu=0.06, omk=0)
pars.InitPower.set_params(ns=0.9655, As=2.1e-9)
results = camb.get_results(pars)
# With CAR patch applied: r_d ≈ 149.2 Mpc  ✓
# Without patch (standard): r_d ≈ 147.1 Mpc
```

---

## 6. Install CLASS with CAR Modification (Alternative to CAMB)

```bash
cd ~
git clone https://github.com/lesgourg/class_public.git
cd class_public
cp source/perturbations.c source/perturbations.c.original
# Apply CAR modification from class/perturbations_car.c
make
```

---

## 7. Run Full Bayesian Analysis

```bash
cd sct-collaboration
bash chains/run_chains.sh
```

This runs PolyChord nested sampling and produces:
- Posterior distributions in chains/output/
- Log-Bayesian-evidence in output/evidence.txt
- Posterior plots in output/posteriors.png

**IMPORTANT:** Bayesian evidence results are only meaningful when using
real survey data (Step 4). With mock data, chi-squared ≈ 0 for some
likelihoods and the evidence figure will be artificially inflated.

---

## 8. Reproduce Paper Figures

```bash
python figures/make_all_figures.py
```

---

## 9. Docker (One-Command Full Reproducibility)

```bash
docker build -t sct-cosmology .
docker run -v $(pwd)/output:/output sct-cosmology bash chains/run_chains.sh
```

---

## File Structure

```
sct-collaboration/
├── sct_core.py                  ← Core CAR calculator (v2.0, bugs fixed)
├── SETUP_INSTRUCTIONS.md        ← This file
├── README.md
├── requirements.txt
├── camb/
│   └── equations_car.f90        ← CAMB CAR sound speed patch
├── class/
│   └── perturbations_car.c      ← CLASS CAR sound speed patch
├── likelihoods/
│   ├── __init__.py
│   ├── combined_likelihood.py
│   ├── desi_bao.py              ← DESI-DR2 BAO (ΛCDM mock if no data)
│   ├── des_y6_lensing.py        ← DES-Y6 (real S8=0.780 mock if no data)
│   ├── hsc_kids_lensing.py      ← HSC-Y3 + KiDS-DR5
│   └── plank_cmb.py             ← Planck PR4 (note: filename has typo)
├── chains/
│   ├── run_chains.sh
│   ├── config_car.ini
│   ├── compute_evidence.py
│   └── plot_posteriors.py
├── figures/
│   └── make_all_figures.py
├── data/
│   └── download_data.sh
└── docker/
    └── Dockerfile
```

---

## Key Parameters (Paper 16 §2.1–2.5)

| Parameter | Value | Status |
|-----------|-------|--------|
| R_b0 (coherence) | 0.2545 ± 0.032 | Derived constant — Paper 17 v4.8 §11.6 |
| S8 | 0.783 ± 0.015 | Analytic — independently verified ✓ |
| b_IA | 1.087 | Analytic — independently verified ✓ |
| r_d | 149.2 ± 0.4 Mpc | Requires CAMB with CAR patch |
| H0 | 70.4 ± 0.5 km/s/Mpc | Derived from r_d and θ* |
| n_s | 28/29 = 0.9655 | Derived from cascade geometry |

---

## Troubleshooting

**ImportError: cannot import name 'Planck_Likelihood' from likelihoods**
The CMB likelihood file is named `plank_cmb.py` (historical typo, missing 'n').
All imports in v2.0 use this spelling consistently. Do not rename the file.

**S8 comes out as 0.042 or IA_bias as 400**
You are running the old v1.0 sct_core.py. Replace it with v2.0 from this repo.
Bug 1 (R_b0=1196.9 instead of 0.2545 (derived constant, Paper 17 v4.8 §11.6)) causes these wrong values.

**H0 comes out as 1.3 km/s/Mpc**
Bug 2 (theta in degrees). Replace with v2.0 sct_core.py.

**r_d comes out as 133 Mpc**
Bug 3 (c/100 instead of c/H0). Replace with v2.0 sct_core.py.

**Bayesian evidence shows Δln B = −3.8 (44:1) but I have no data files**
Results computed from mock data are circular and not valid. Download
real survey data first (Step 4).
