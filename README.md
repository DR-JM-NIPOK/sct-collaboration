# class/ — Modified CLASS Boltzmann Solver for CAR

## Overview

This directory contains the patch to implement the Codified Acoustic Relation
(CAR) sound speed within the CLASS (Cosmic Linear Anisotropy Solving System)
Boltzmann solver.

**The single modification to CLASS:**

| | Formula |
|---|---|
| Standard ΛCDM | `cs2 = 1.0 / (3.0 * (1.0 + R))` |
| **CAR** | **`cs2 = (1.0 + R) / 3.0`** |

CLASS is an alternative to CAMB — both produce identical CAR results when
the sound speed is patched correctly. CLASS is preferred for:
- Faster execution on certain architectures
- Built-in Boltzmann hierarchy output
- Integration with MontePython

## Installation

```bash
# 1. Clone CLASS
git clone https://github.com/lesgourg/class_public.git
cd class_public

# 2. Apply the CAR patch
patch -p1 < ../sct-collaboration/class/perturbations_CAR.patch

# 3. Build
make clean && make -j4

# 4. Install Python wrapper
cd python
python setup.py install

# 5. Verify
python ../sct-collaboration/class/class_car_test.py
```

## Expected test output

```
CLASS-CAR verification
  Standard CLASS r_s : 150.00 Mpc
  CAR-modified r_s   : 149.09 Mpc
  Δr_s               : -0.91 Mpc
  c_s²(z_dec)        : 0.4107  (CAR)  vs  0.2792  (ΛCDM)
  PASS
```

## CLASS version compatibility

Tested with CLASS v2.9.x and v3.x. The sound speed is computed in
`source/perturbations.c` in function `perturb_derivs()`.

## Files

| File | Description |
|---|---|
| `perturbations_CAR.patch` | Unified diff for `perturbations.c` |
| `class_car_test.py` | Verification test |
| `README.md` | This file |
