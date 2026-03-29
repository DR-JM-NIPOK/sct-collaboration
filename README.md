# figures/ — Paper Figure Reproduction Scripts

## Overview

This directory contains scripts to reproduce all figures from
SCT Cosmology Series Paper #16.

## Quick Start

```bash
# Generate all figures (Figs 1, 2, 5 require no chain output)
python figures/make_all_figures.py --output output/figures/

# Single figure
python figures/make_all_figures.py --fig 1

# High-resolution for publication
python figures/make_all_figures.py --format pdf --dpi 300
```

## Figure Index

| Fig | File | Description | Requires |
|---|---|---|---|
| 1 | `fig1_tension_summary.py` | H₀, S₈, r_d tension summary | `sct_core.py` only |
| 2 | `fig2_sound_horizon.py` | c_s(z) and r_d integrand | `sct_core.py` only |
| 3 | `fig3_posteriors.py` | 2D posterior contours | PolyChord chains |
| 4 | `fig4_residuals.py` | Standardised residuals | Likelihood modules |
| 5 | `fig5_forecasts.py` | Euclid/LSST/CMB-S4 kill switch | `sct_core.py` only |

Figures 1, 2, and 5 can be generated immediately with:
```bash
pip install numpy scipy matplotlib
python figures/make_all_figures.py
```

Figures 3 and 4 require completed PolyChord chains and data.
See `chains/README.md` for instructions.

## Expected Output Files

```
output/figures/
├── fig1_tension_summary.pdf
├── fig2_sound_horizon.pdf
├── fig3_posteriors.pdf         (requires chains)
├── fig4_residuals.pdf          (requires data + chains)
└── fig5_forecasts.pdf
```
