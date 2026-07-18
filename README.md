# sct-collaboration v4.8.5

**Successive Collision Theory (SCT) — public reference implementation**

DR JM NIPOK · N.J.I.T. · ORCID [0009-0006-3940-4450](https://orcid.org/0009-0006-3940-4450)
Paper 15 DOI: [10.13140/RG.2.2.10321.29288](https://doi.org/10.13140/RG.2.2.10321.29288)
Series 2 Paper 1 DOI: [10.13140/RG.2.2.14355.03366](https://doi.org/10.13140/RG.2.2.14355.03366)
License: GPL-3.0

---

## What is in this release

This is the post-NLA-recursive-audit release of the SCT computational
framework. Every file has been brought into consistency with the canonical
Series 2 Paper 1 derivation of R_b, and the audit findings are documented
in [CHANGELOG.md](CHANGELOG.md).

### Status snapshot (verified by audit, April 2026)

| Quantity | Value | Status |
|---|---|---|
| R_b (derived, §11.6) | 0.2545 ± 0.032 | Derived from SO(3) cascade + QCD junction |
| c_s² (z=0) | 0.41817 ± 0.011 | Derived from (1 + R_b)/3 |
| b_IA | 1.0848 ± 0.011 | Matches IA analyses (DES-Y6, KiDS-DR5, HSC-Y3) |
| S8 (analytic) | 0.7988 | From 0.832 × (1 + R_b/3)^(-1/2) |
| S8 (numerical) | 0.7838 ± 0.015 | Matches DES-Y6 0.780 ± 0.012 within 0.4σ |
| r_d (standard horizon) | 146.8 ± 5 Mpc | Consistent with DESI-DR2 BAO (147 ± 1); CAR does not modify it |
| N_eff (predicted) | 2.514 ± 0.05 | CMB-S4 testable at 17.7σ vs SM N_eff = 3.046 [^neff] |

### r_d: standard recombination horizon (RNLA v2.3)

The CAR enhanced sound speed c_s²(z) = (1 + R_b(z))/3 is a **late-time**
coherent-sector quantity (it sets S8 and b_IA); it does **not** modify the
recombination acoustics. The recombination/drag sound horizon is therefore
**standard**, r_d ≈ 146.8 Mpc (r_*(z*) ≈ 144.4 Mpc), **consistent with
DESI-DR2 BAO (147 ± 1 Mpc)**. The v4.8.1 value 161.4 Mpc was a category
error (late-time c_s applied at recombination) and is retired. The S8 and b_IA predictions
remain robust and observationally consistent. Paper 15 (v2.22, 2026-06) already
reflects this correction. See `CHANGELOG.md` for full audit details.

[^neff]: **Open register question (2026-07-10 in-session verification):** naively
substituting N_eff = 2.514 into stock ΛCDM at fixed θ* is excluded at Δχ² ≈ +658
against Planck TT (83 binned points); the register's H0≈66-67 note reproduces only
under standard N_eff=3.044, not under N_eff=2.514. This N_eff prediction requires
the CAR-patched Boltzmann run to confirm that SCT's 2.514 and the quantity CMB-S4
will actually measure are the same observable. See `tensions.csv` T031 and
`predictions.csv` P088.

---

## Repository structure

```
sct-collaboration/
├── sct_core.py              ← Canonical CAR core calculator (v4.8.4)
├── predictions.csv          ← All SCT predictions w/ falsification criteria
├── tensions.csv             ← Observational tensions reference
├── bh_interior.py           ← Black hole interior structure (Paper 16)
├── coherence.py             ← Coherence cascade calculations (Paper 5)
├── growth.py                ← Structure growth (Paper 8)
├── hereditary.py            ← Hereditary mode propagation (Paper 12)
├── camb/                    ← CAMB modification (corrected v4.8.4)
│   ├── equations_car.f90
│   ├── equations_CAR.patch
│   └── equations_car_test.py
├── class/                   ← CLASS modification (corrected v4.8.4)
│   ├── perturbations_CAR.patch
│   └── class_car_test.py
├── chains/                  ← MCMC / nested-sampling configurations
├── likelihoods/             ← Per-survey likelihoods (DESI, DES, HSC, KiDS, Planck)
├── data/                    ← Mock data generator
├── figures/                 ← Figure generation scripts
├── docker/                  ← Reproducible Docker environment
├── tests/                   ← Pytest suite
└── paper17-v4-deriver-Rb/   ← Series 2 Paper 1 §11.6 R_b derivation reference
```

---

## Quick start

### Verify canonical predictions (no CAMB required)

```bash
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration
pip install -r requirements.txt
python3 sct_core.py                # compact summary
python3 sct_core.py --validate     # full validation report with σ values
python3 sct_core.py --verbose      # detailed intermediate values
```

Expected output highlights:
```
R_b0  (DERIVED §11.6)      0.2545
c_s² = (1+R_b)/3           0.41817
b_IA  = 1 + R_b/3          1.08483
S8 (analytic)              0.7988
S8 (numerical, −0.015)     0.7838
r_d  (standard, Mpc)       146.80       [consistent with DESI-DR2 BAO]
N_eff (SCT predicted)      2.514        CMB-S4 separation 17.7σ
```

### Apply the CAMB CAR patch (optional — for full Boltzmann run)

```bash
cd <your-camb-source-dir>
patch -p1 < /path/to/sct-collaboration/camb/equations_CAR.patch
make                                    # rebuild CAMB
cd /path/to/sct-collaboration
python3 camb/equations_car_test.py     # verify
```

The recombination sound horizon is **standard** (r_drag ≈ 146.8 Mpc) and is
**not** modified by SCT — do not patch CAMB's recombination cs2. The CAR
coherent sound speed is a late-time effect (S8, b_IA); see
`camb/equations_car.f90` (module `CAR_corrections`). The test script confirms
the standard r_drag ≈ 146.8 Mpc.

### Run the test suite

```bash
pytest tests/ -v
```

### Docker (reproducible environment)

```bash
docker build -t sct-coll:v4.8.4 -f docker/Dockerfile.txt .
docker run --rm -v $(pwd):/workspace sct-coll:v4.8.4 python3 sct_core.py --validate
```

---

## Papers in the SCT Cosmology Series

The repository implements computational tools for Paper 1–46 of the
SCT Cosmology Series. Direct dependencies:

- **Paper 15** — From Chaos to Codified Acoustics (CAR ansatz; r_d, S8, H0 predictions)
- **Series 2 Paper 1** — From Chaos to Closure (Section 11.6 derives R_b = 0.2545)
- **Paper 9** — From Chaos to Convergence (QCD junction conditions, 13.6% loss)
- **Paper 5** — From Chaos to Cascading Coherence (SO(3) angular momentum)

Full series: https://thenaturalstateofnature.org/PREPRINTS/From_Chaos_To_Consilience/

---

## v4.8.1 audit acknowledgment

The NLA recursive audit that produced this release identified eleven
inconsistencies between canonical Python (sct_core.py at root) and
the rest of the repository — including two critical bugs in the CAMB
and CLASS Fortran/C patches that produced unphysical (superluminal)
sound speeds at z=0. All findings are documented in `CHANGELOG.md`.
The audit did not change any derived theoretical constants from
Series 2 Paper 1 §11.6; it only ensured every file in the repository
implements them correctly and reports outputs that the code can
actually reproduce.

---

## License & citation

GPL-3.0. If you use this code in published work, please cite:

```bibtex
@misc{nipok2026sct,
  author       = {{DR JM NIPOK}},
  title        = {SCT Cosmology Series (Paper 1--46) Computational Framework},
  year         = {2026},
  publisher    = {GitHub},
  howpublished = {\url{https://github.com/DR-JM-NIPOK/sct-collaboration}},
  doi          = {10.13140/RG.2.2.10321.29288}
}
```
