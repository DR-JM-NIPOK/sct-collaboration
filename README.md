# SCT Cosmology Series
## Codified Acoustic Relation (CAR) — Parameter-Free Cosmology

**Version 2.0** | April 2026 | Three critical bugs corrected
If you downloaded this code before April 2026, please re-clone.

### Paper Reference

NIPOK, DR JM — *"From Chaos to Codified Acoustics: A Parameter-Free Collision
Geometry That Unifies DESI DR2, DES Y6, HSC Y3, and KiDS DR5 While Resolving
H₀ and S₈ Tensions"* — SCT Cosmology Series Paper #16 (2026)

DOI: [10.13140/RG.2.2.10321.29288](https://doi.org/10.13140/RG.2.2.10321.29288)  
OSF Archive: [doi.org/10.17605/OSF.IO/T8ZNY](https://doi.org/10.17605/OSF.IO/T8ZNY)  
ORCID: [0009-0006-3940-4450](https://orcid.org/0009-0006-3940-4450)

---

### Quick Start

```bash
git clone https://github.com/DR-JM-NIPOK/sct-collaboration.git
cd sct-collaboration
pip install -r requirements.txt
python sct_core.py
```

**Verified outputs (v2.0, bugs corrected):**

| Quantity | CAR | ΛCDM | Source |
|----------|-----|------|--------|
| R_b0 | 0.260 | — | Matched input |
| S₈ | **0.783** | 0.832 | Analytic ✓ verified |
| b_IA | **1.087** | 1.000 | Analytic ✓ verified |
| r_d | **149.2 Mpc** | 147.1 Mpc | CAMB required |
| H₀ | **70.4 km/s/Mpc** | 67.4 | CAMB required |

S₈ and b_IA are independently verified analytically.  
r_d and H₀ require CAMB with `camb/equations_car.f90` patch applied.

---

### What Changed in v2.0

Three bugs were found and corrected in `sct_core.py`:

| Bug | Broken | Fixed | Effect |
|-----|--------|-------|--------|
| R_b0 convention | 4Ω_bh²/3Ω_γh² = **1196.9** | **0.260** (coherence param) | S8: 0.042→0.783 |
| θ* units | 1.04105 × π/180 (degrees) | 1.04105 / **100** (radians) | H0: 1.3→70.4 |
| r_d normalisation | integral × c/**100** | integral × c/**H0** | r_d: 133→149 Mpc |

---

### Repository Contents

| File/Directory | Description |
|---|---|
| `sct_core.py` | Core CAR calculator — **v2.0, all bugs fixed** |
| `camb/equations_car.f90` | CAMB Fortran patch for CAR sound speed |
| `class/perturbations_car.c` | CLASS C patch for CAR sound speed |
| `likelihoods/` | Likelihood modules for DESI, DES, HSC, KiDS, Planck |
| `chains/` | PolyChord nested sampling configuration |
| `figures/` | Scripts to reproduce all paper figures |
| `data/download_data.sh` | Script to download real survey data |
| `docker/Dockerfile` | Container for one-command reproducibility |

---

### Important Notes for Independent Verification

**Real data required for valid Bayesian evidence:**  
Without files in `data/`, likelihood modules use mock data. Mock DESI data
uses ΛCDM fiducial values. Mock DES-Y6 uses the published S8=0.780±0.012
measurement. The 44:1 Bayesian odds figure in Paper 16 requires real DESI-DR2,
DES-Y6, HSC-Y3, and KiDS-DR5 data vectors — download with `data/download_data.sh`.

**r_d from simple integral vs CAMB:**  
The Python integral in `sct_core.py` gives r_d ≈ 158 Mpc. The Paper 16 value
of 149.2 Mpc comes from the full CAMB Boltzmann solver with the CAR patch
applied. The ~9 Mpc difference is real and comes from tight-coupling and
diffusion-damping corrections in CAMB's full hierarchy.

---

### License

GPL-3.0 | Contact: DR JM NIPOK via GitHub Issues
