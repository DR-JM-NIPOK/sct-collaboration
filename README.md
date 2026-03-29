# data/ — Data Vectors and Download Instructions

## Overview

This directory contains scripts for obtaining the observational data used
in SCT Cosmology Series Paper #16. The data themselves are not distributed
here because they are the property of the respective survey collaborations
and are available from their official data portals.

## Data Sources

| Survey | Data Type | n_data | URL | Reference |
|---|---|---|---|---|
| DESI-DR2 | BAO measurements | 32 | https://data.desi.lbl.gov/public/dr2/bao/ | arXiv:2603.12458 |
| DES-Y6 | 3×2pt correlation functions | 1800 | https://des.ncsa.illinois.edu/releases/y6 | arXiv:2603.08941 |
| HSC-Y3 | Cosmic shear ξ±(θ) | 576 | https://hsc-release.mtk.nao.ac.jp | ApJ 960, 62 (2025) |
| KiDS-DR5 | Cosmic shear | 912 | https://kids.strw.leidenuniv.nl/DR5/ | A&A 688, A120 (2026) |
| Planck PR4 | CMB TT/TE/EE | 86 (compressed) | https://pla.esac.esa.int | A&A 687, A160 (2024) |

## Mock Data (for testing)

To test the pipeline without real data:

```bash
python data/mock_data_generator.py --output data/ --noise 0.02
```

This generates HDF5 mock files at the CAR best-fit point:

```
data/
├── desi_dr2_bao.hdf5      (32 mock DV/r_d measurements)
├── des_y6_3x2pt.hdf5      (1800 mock 3×2pt data points)
├── hsc_y3_wl.hdf5         (576 mock shear measurements)
├── kids_dr5_wl.hdf5       (912 mock shear measurements)
└── planck_pr4_compressed.hdf5  (compressed CMB likelihood)
```

## Downloading Real Data

```bash
# Print download instructions
python data/mock_data_generator.py --real-data
```

## HDF5 File Format

All likelihood modules expect HDF5 files with the following structure:

**DESI-DR2:**
```
/data_vector    shape: (32,)   — DV(z)/r_d values
/covariance     shape: (32,32) — full covariance matrix
attrs: survey, n_bins, z_eff_values
```

**DES-Y6 / HSC / KiDS:**
```
/xi_plus        shape: (n_xi,)   — ξ+(θ) shear correlation
/xi_minus       shape: (n_xi,)   — ξ-(θ) shear correlation
/gamma_t        shape: (n_ggl,)  — γ_t galaxy-galaxy lensing
/w_theta        shape: (n_clust,)— w(θ) clustering (DES only)
/covariance     shape: (n,n)     — full joint covariance
attrs: survey, n_z_bins, theta_min_arcmin, theta_max_arcmin
```

**Planck PR4 (compressed):**
```
/compressed_params   shape: (86,)   — compressed CMB parameters
/fisher_matrix       shape: (86,86) — inverse covariance (Fisher)
attrs: lmax, likelihood_version
```
