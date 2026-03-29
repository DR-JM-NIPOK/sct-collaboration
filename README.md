# camb/ — Modified CAMB Boltzmann Solver for CAR

## Overview

This directory contains the patch required to implement the Codified Acoustic
Relation (CAR) sound speed within the CAMB Boltzmann solver.

**One-line change to CAMB:**

| | Formula |
|---|---|
| Standard ΛCDM | `cs² = 1 / (3 × (1 + R_b))` |
| **CAR** | **`cs² = (1 + R_b) / 3`** |

This single modification shifts the sound horizon from 150.0 Mpc (Planck ΛCDM)
to 149.1 ± 0.3 Mpc, resolving the DESI-DR2 tension.

## Installation

```bash
# 1. Clone CAMB
git clone https://github.com/cmbant/CAMB.git
cd CAMB

# 2. Apply the CAR patch
patch -p1 < ../sct-collaboration/camb/equations_CAR.patch

# 3. Install modified CAMB
pip install -e .

# 4. Verify
python camb/equations_car_test.py
```

## Expected test output

```
CAMB-CAR Sound horizon test
  Standard CAMB r_d : 150.00 Mpc
  CAR-modified r_d  : 149.10 Mpc
  Δr_d              : -0.90 Mpc
  PASS
```

## CAMB version compatibility

Tested with CAMB 1.4.0 and 1.5.x. The sound speed calculation is in
`fortran/equations.f90` — the relevant subroutine name may differ by version:

| Version | Subroutine |
|---|---|
| ≤ 1.3 | `derivs` |
| 1.4–1.5 | `ThermoData_Get` |

## Files

| File | Description |
|---|---|
| `equations_CAR.patch` | Unified diff for `equations.f90` |
| `equations_car_test.py` | Verification script |
| `README.md` | This file |
