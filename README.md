# chains/ — Nested Sampling Configuration and Analysis

## Overview

This directory contains all configuration and analysis scripts for running
PolyChord nested sampling to:
1. Compute posteriors for CAR and ΛCDM parameters
2. Compute Bayesian evidence for model comparison
3. Generate posterior plots using GetDist

## Directory Contents

| File | Description |
|---|---|
| `run_chains.sh` | Master bash script — start here |
| `run_polychord.py` | PolyChord Python runner with prior transforms |
| `config_car.ini` | CAR model settings (2 free parameters) |
| `config_lcdm.ini` | ΛCDM model settings (6 + 42 nuisance parameters) |
| `compute_evidence.py` | Extract log-evidence and compute Bayes factors |
| `plot_posteriors.py` | GetDist posterior contour plots |
| `README.md` | This file |

## Prerequisites

```bash
# PolyChord (requires MPI + Fortran compiler)
git clone https://github.com/PolyChord/PolyChordLite.git
cd PolyChordLite
make
pip install -e .

# GetDist (for posterior plots)
pip install getdist
```

## Quick Start

```bash
# Step 1: Generate mock data (if real data not available)
python data/mock_data_generator.py --output data/

# Step 2: Run CAR nested sampling
bash chains/run_chains.sh --model car --data combined

# Step 3: Run ΛCDM for comparison
bash chains/run_chains.sh --model lcdm --data combined

# Step 4: Compute Bayes factor
python chains/compute_evidence.py \
    --compare-car  output/evidence_car_combined.txt \
    --compare-lcdm output/evidence_lcdm_combined.txt \
    --output       output/bayes_factor.txt
```

## Expected Outputs

After a full combined run, the key numbers to verify:

| Quantity | Expected Value |
|---|---|
| ln Z (CAR) | −1248.52 ± 0.04 |
| ln Z (ΛCDM) | −1247.32 ± 0.03 |
| Δln B (DESI+Planck) | −1.20 ± 0.28 (5:1 odds) |
| Δln B (conservative combined) | −3.80 ± 0.40 (44:1 odds) |
| CAR Ω_m posterior | 0.312 ± 0.009 |
| CAR R_b (derived) | 0.257 ± 0.032 |  (Paper 17 v4.0 Section 11.6 — not sampled)

## Model Comparison Notes

The Bayesian Information Criterion (BIC) comparison uses:

```
BIC_ΛCDM = χ²_ΛCDM + k_ΛCDM × ln(N)   k_ΛCDM = 48 (6 cosmo + 42 nuisance)
BIC_CAR  = χ²_CAR  + k_CAR  × ln(N)    k_CAR  = 2
ΔBIC = BIC_CAR - BIC_ΛCDM ≈ −346.6     (decisive evidence for CAR)
```

The conservative log-evidence from Bayesian dimensionality correction:
```
Δln B_conservative = −3.80 ± 0.40    (odds 44:1)
```

See §8 of the paper and `compute_evidence.py` for full derivation.

## Runtime Estimates

| Configuration | CPUs | Estimated Time |
|---|---|---|
| CAR (2 params, n_live=500) | 1 | ~2 hours |
| ΛCDM (48 params, n_live=500) | 16 | ~48 hours |
| Full combined (CAR) | 4 | ~8 hours |

Use Docker for a reproducible environment with all dependencies.
