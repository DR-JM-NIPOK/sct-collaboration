# camb/ — Modified CAMB Boltzmann Solver for CAR

## Overview

This directory contains the patch required to implement the Codified Acoustic
Relation (CAR) sound speed within the CAMB Boltzmann solver.

**One-line change to CAMB:**

| | Formula |
|---|---|
| Standard ΛCDM | `cs² = 1 / (3 × (1 + R_b))` |
| **CAR** | **`cs² = (1 + R_b) / 3`** |

This modification is a **late-time** coherent-sector effect (it sets S8 and
b_IA). It does **not** modify the recombination sound horizon, which stays
**standard** (r_drag ≈ 146.8 Mpc, consistent with DESI-DR2 BAO 147 ± 1 Mpc).
RNLA v2.3 (June 2026): the earlier "shifts r_d to 149.1/161.4 Mpc" claim was a
category error (late-time c_s applied at recombination) and is retired.

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
  Standard CAMB r_drag : 146.80 Mpc   (recombination horizon — UNCHANGED by CAR)
  r_*(z*)              : 144.40 Mpc
  DESI-DR2 r_d         : 147.0 ± 1.0 Mpc   (consistent, ~0.2 sigma)
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
