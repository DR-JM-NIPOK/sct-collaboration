# SCT Cosmology — Complete Setup Instructions

## Overview

This document describes how to reproduce the full analysis from
SCT Cosmology Series Paper #16 (CAR framework).

**Three levels of setup:**

| Level | Requirements | Time | What you get |
|---|---|---|---|
| **Basic** | Python + NumPy/SciPy | 2 min | Core CAR predictions |
| **Analysis** | + CAMB or CLASS | 30 min | Full Boltzmann solver |
| **Full** | + PolyChord + MPI | 1-2 days | Posteriors + Bayes evidence |

---

## Level 1: Basic (Core CAR Calculator)

```bash
# Clone the repository
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration

# Install minimum dependencies
pip install numpy scipy matplotlib

# Run the CAR calculator
python sct_core.py
```

**Expected output:**
```
========================================================
  Codified Acoustic Relation (CAR) — Predictions
========================================================
  Quantity               CAR          ΛCDM         Δ
--------------------------------------------------------
  R_b0                  0.2600       0.2600    (same)
  c_s²(z→∞)  [×c²]     0.42000      0.27895  +0.14105
  r_d  [Mpc]           149.10       150.00     -0.90
  H₀  [km/s/Mpc]        70.40        67.40     +3.00
  S₈  (numerical)        0.783        0.832    -0.049
  b_IA                   1.087        1.000    +0.087
```

---

## Level 2: Full Python Analysis

```bash
# Install all Python dependencies
pip install -r requirements.txt

# Generate mock data (real data requires survey registration)
python data/mock_data_generator.py --output data/

# Run unit tests
pytest tests/ -v

# Generate Figures 1, 2, and 5 (no chains needed)
python figures/make_all_figures.py --output output/figures/
```

---

## Level 3: Boltzmann Solvers

### CAMB (recommended)

```bash
# Install unmodified CAMB first
pip install camb

# Apply the CAR sound speed patch
cd /path/to/CAMB_source   # if building from source
patch -p1 < /path/to/sct-collaboration/camb/equations_CAR.patch
make clean && make

# Verify
python sct-collaboration/camb/equations_car_test.py
```

### CLASS (alternative)

```bash
# Clone CLASS
git clone https://github.com/lesgourg/class_public.git
cd class_public

# Apply CAR patch
patch -p1 < /path/to/sct-collaboration/class/perturbations_CAR.patch

# Build
make clean && make -j4
cd python && python setup.py install

# Verify
python /path/to/sct-collaboration/class/class_car_test.py
```

---

## Level 4: Nested Sampling (Full Posterior + Evidence)

### Prerequisites

```bash
# System dependencies (Ubuntu/Debian)
sudo apt install -y gfortran gcc g++ make cmake openmpi-bin libopenmpi-dev

# PolyChord (nested sampling)
git clone https://github.com/PolyChord/PolyChordLite.git
cd PolyChordLite
make
pip install -e .

# GetDist (posterior plots)
pip install getdist
```

### Running the Analysis

```bash
# Generate mock data (or download real data — see data/README.md)
python data/mock_data_generator.py --output data/

# Run CAR nested sampling (full combined: ~8 hrs on 4 cores)
bash chains/run_chains.sh --model car --data combined

# Run ΛCDM for comparison (~48 hrs on 16 cores)
bash chains/run_chains.sh --model lcdm --data combined

# Compute Bayes factor
python chains/compute_evidence.py \
    --compare-car  output/evidence_car_combined.txt \
    --compare-lcdm output/evidence_lcdm_combined.txt \
    --output output/bayes_factor_combined.txt

# Generate all figures
python figures/make_all_figures.py --output output/figures/
```

---

## Docker (Complete Reproducibility)

```bash
# Build container (includes all dependencies except PolyChord)
docker build -t sct-cosmology docker/

# Run CAR calculator
docker run sct-cosmology

# Run figures
docker run -v $(pwd)/output:/output sct-cosmology \
    python figures/make_all_figures.py --output /output/figures/

# Full analysis (requires PolyChord in container — see docker/Dockerfile)
docker-compose -f docker/docker-compose.yml run sct-car
```

---

## Expected Results Summary

After running the full analysis:

| Quantity | CAR Prediction | Observation | Tension |
|---|---|---|---|
| r_d | 149.1 ± 0.3 Mpc | DESI-DR2: 147.0 ± 1.0 | 2.1σ |
| H₀ | 70.4 ± 0.4 km/s/Mpc | SH0ES: 73.0 ± 1.0 | 2.4σ |
| S₈ | 0.783 ± 0.015 | DES-Y6: 0.780 ± 0.012 | 0.2σ |
| b_IA | 1.087 ± 0.002 | DES-Y6: 1.08 ± 0.04 | 0.2σ |
| ln Z (CAR) | −1248.52 ± 0.04 | — | — |
| Δln B | −3.80 ± 0.40 | — | 44:1 odds for CAR |

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'sct_core'`**
→ Run from the `sct-collaboration/` root directory, or:
```bash
export PYTHONPATH=/path/to/sct-collaboration:$PYTHONPATH
```

**PolyChord build fails**
→ Ensure `gfortran` is installed and MPI is configured.
```bash
which mpif90 && mpif90 --version
```

**CAMB returns r_d ≈ 150 Mpc (not 149.1)**
→ The CAR patch was not applied. Re-apply `equations_CAR.patch` and rebuild.

**Tests fail on `test_conservative_evidence`**
→ Requires `chains/compute_evidence.py` to be importable. Run from repo root.

---

## File Structure After Full Setup

```
sct-collaboration/
├── sct_core.py              ← Core CAR calculator ✓
├── requirements.txt         ✓
├── camb/                    ← Patched CAMB (Level 3)
├── class/                   ← Patched CLASS (Level 3)
├── likelihoods/             ← Survey likelihood modules ✓
├── chains/                  ← PolyChord analysis (Level 4)
├── figures/                 ← Figure scripts ✓
├── data/                    ← Mock data generator ✓
├── docker/                  ← Container definition ✓
├── tests/                   ← Unit tests ✓
└── output/                  ← Generated during analysis
    ├── chains_car_combined/
    ├── chains_lcdm_combined/
    ├── figures/
    └── evidence_*.txt
```
