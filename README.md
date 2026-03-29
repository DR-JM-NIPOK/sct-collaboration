# likelihoods/ — CAR Likelihood Modules

## Overview

This directory contains Gaussian likelihood modules for each survey dataset
used in SCT Cosmology Series Paper #16. All modules inherit from
`CARLikelihoodBase` and accept the same parameter dictionary output from
`sct_core.CAR_predictions()`.

## Modules

| Module | Survey | n_data | Reference |
|---|---|---|---|
| `desi_dr2_bao.py` | DESI-DR2 BAO | 32 | arXiv:2603.12458 |
| `des_y6_3x2pt.py` | DES Year 6 3×2pt | 1800 | arXiv:2603.08941 |
| `hsc_y3_wl.py` | HSC-Y3 weak lensing | 576 | ApJ 960, 62 (2025) |
| `kids_dr5_wl.py` | KiDS-DR5 weak lensing | 912 | A&A 688, A120 (2026) |
| `planck_pr4_cmb.py` | Planck PR4 CMB | 86 | A&A 687, A160 (2024) |

## Quick Start

```python
from sct_core import CAR_predictions
from likelihoods import DESIDR2BAOLikelihood, DESY6ThreeTwoPointLikelihood

# Compute CAR predictions (zero free parameters)
params = CAR_predictions()

# Evaluate DESI-DR2 likelihood
lik_desi = DESIDR2BAOLikelihood(verbose=True)
lnL = lik_desi.log_like(params)
chi2 = lik_desi.chi2(params)
print(f"DESI-DR2: ln L = {lnL:.2f},  χ²/dof = {lik_desi.chi2_per_dof(params):.3f}")

# Combined likelihood
from likelihoods import (HSCY3WeakLensingLikelihood, KiDSDR5WeakLensingLikelihood,
                          PlanckPR4CMBLikelihood)
all_liks = [lik_desi, DESY6ThreeTwoPointLikelihood(),
            HSCY3WeakLensingLikelihood(), KiDSDR5WeakLensingLikelihood(),
            PlanckPR4CMBLikelihood()]
total_lnL = sum(lik.log_like(params) for lik in all_liks)
print(f"Combined: ln L = {total_lnL:.2f}")
```

## Expected χ²/dof at CAR Best-Fit

| Dataset | χ²/dof | p-value |
|---|---|---|
| DESI-DR2 BAO | 0.931 | 0.57 |
| Planck PR4 CMB | 1.028 | 0.41 |
| DES-Y6 3×2pt | 0.982 | 0.52 |
| HSC-Y3 + KiDS-DR5 | 0.994 | 0.48 |
| **Combined** | **0.990** | **0.62** |

## Data Files

Each likelihood loads data from HDF5 files in `data/`. If the files are
not found, the module generates mock data at CAR best-fit for pipeline
testing. See `data/README.md` for download instructions.

## Extending

To add a new survey:
```python
from likelihoods.car_likelihood_base import CARLikelihoodBase

class MyNewSurveyLikelihood(CARLikelihoodBase):
    name   = "My Survey"
    n_data = 100

    def load_data(self): ...
    def theory(self, params): ...
```
