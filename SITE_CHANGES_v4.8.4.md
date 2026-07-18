# sct-collaboration v4.8.4 — change record (RNLA v2.3 recursive audit)

Audit: Recursive NIPOK Level Audit (RNLA) v2.3 — standing checks 13 (cross-observable
consistency), 14 (repository <-> paper consistency), 15 (post-diction discipline).
Base: v4.8.3. All branches audited (root, camb, chains, class, data, docker, figures,
likelihoods, paper17-v4-deriver-Rb, tests).

## Primary fix — r_d category error (checks 13/14)
The v4.8.1 audit had set r_d = 161.4 Mpc and declared the earlier 146.8 Mpc "retracted."
RNLA v2.3 reverses this: **161.4 Mpc was itself the category error.**

- The CAR enhanced sound speed cs2 = (1 + R_b)/3 (R_b0 = 0.2545) is a **late-time
  coherent-sector** quantity — it sets S8 and the IA bias b_IA. It is **not** the
  recombination acoustic speed. At recombination R_b -> 0, so (1 + R_b)/3 -> 1/3, the
  pure-radiation speed; omitting baryon loading inflates the horizon to 161.4/165/178 Mpc.
- The recombination/drag sound horizon is the **standard photon-baryon horizon**, fixed
  by omega_b and omega_m (H0-independent): **r_d (drag) = 146.8 +/- 5 Mpc, r_*(z*) = 144.4 Mpc.**
- Numerically verified (standard baryon-loaded speed): r_drag = 146.97 Mpc, r_* = 144.31 Mpc
  at Planck cosmology — matching the canonical 146.8/144.4 and Planck 2018.

## Consequences
- **DESI-DR2 BAO**: standard r_d = 146.8 Mpc is **consistent** with DESI (147 +/- 1 Mpc, ~0.2 sigma).
  SCT introduces no early-time BAO tension. (Prior "~14 sigma tension / does not close DESI" retired.)
- **H0**: global theta* + r_d = 146.8 gives **H0 ~ 66.3 km/s/Mpc** (PARTIAL). The prior
  70.4 km/s/Mpc is **not CMB-derivable** and is retired; local H0 is raised toward 70-73 by
  the late-time void + temporal mechanism. (Matches RNLA check 13: compute_H0_from_theta_and_rd -> 66.3.)

## Files changed (branch copies kept byte-identical where shared)
- sct_core.py x10: R_D_DERIVED 161.4->146.8, R_D_UNCERTAINTY 0.3->5.0, added R_STAR_DERIVED=144.4;
  compute_r_d_integral now uses the standard baryon-loaded speed; CAR cs scoped late-time;
  H0 self-consistent from theta*+r_d (66.3); headers/reports/CLI rewritten; version -> v4.8.4.
- camb/equations_car.f90, class/perturbations_car.c: reframed to **late-time** scope; recombination
  cs2 left standard; the doubly-wrong density-ratio-in-CAR-formula form removed from the C file.
- camb/equations_CAR.patch, class/perturbations_CAR.patch: **retired** (deprecation notices, no hunks).
- camb/equations_car_test.py, class/class_car_test.py: expectations -> standard 146.8.
- tests/test_sct_core.py, tests/test_likelihoods.py: r_d / tension / monotonicity tests corrected;
  suite passes 54/54.
- audit_framework.py: verify_r_d_chain now verifies the standard 146.8 horizon and isolates the
  161.4 category error; report title -> v4.8.4.
- likelihoods/desi_dr2_bao.py, desi_bao.py, combined_likelihood.py: r_d narrative -> 146.8 / consistent.
- tensions.csv x2: T006 (r_d) -> consistent 146.8; T001 (H0) -> 70.4 retired, global ~66.3; T003 (A_lens)
  -> accommodation/post-diction (check 15).
- SETUP_INSTRUCTIONS.md x9, READMEs (root, camb, docker, paper17): stale example outputs corrected
  (R_b=0.2545, c_s2=0.41817, b_IA=1.0848, r_d=146.8, H0=66.3).
- CHANGELOG.md x10: [4.8.4] entry appended (append-only history).

## Module E — bh_interior.py (TOV + QCD floor)
- Replaced a non-convergent brentq EOS inversion (bracket 1e10-1e36, xtol too tight) with the
  exact analytic inverse of the piecewise EOS (machine precision).
- Fixed two dimensional bugs: a stray c^2 in the EOS normalization (P was ~c^2*eps, superluminal)
  and in the dP/dr prefactor (pressure scale height was ~1e-10 m). TOV is now the standard SI form
  dP/dr = -(G/c^4)(eps+P)(Mc^2+4pi r^3 P)/[r^2(1-2GM/c^2 r)]; dM/dr = 4pi r^2 eps/c^2.
- Module now runs end-to-end and yields km-scale neutron stars. With the current stiff EOS
  normalization M_max ~ 3.0-3.3 M_sun at R~18 km; exact M_max/R is an EOS-calibration choice
  (P/eps at saturation, Gamma, cs^2) flagged in the file docstring.

## A_lens (check 15)
predictions.csv P007 was already correct (post-diction, CONSISTENT, ~1.18). tensions.csv T003
brought into agreement: A_lens is an **accommodation of the pre-existing Planck anomaly**, labelled
post-diction — never "derives 1.19" / "resolved."

## Unchanged (verified still canonical)
A* = 6.173, N_coh = 14.06, f_b = 0.162 (6.17 re-cascade); R_b = 0.2545, c_s2 = 0.41817,
b_IA = 1.0848, N_eff = 2.514; predictions catalog = 115; CAMB/CLASS cosmology sector
(R_b, c_s2, N_eff) untouched. No residual 5.970 / 13.51 / 0.1675 in live code.

## Verification
- python3 -m pytest tests/ -> 54 passed.
- python3 audit_framework.py -> 5 PASS, 0 FAIL (standard r_d = 147.0 Mpc within canonical 146.8 +/- 5).
- Stale-value sweep (R_b 0.257, b_IA 1.087, cs2 0.419/0.420, A* 5.970, N_coh 13.51, f_b 0.1675,
  r_d 161.4/149.1, H0 70.4) -> no live occurrences; remaining mentions are dated CHANGELOG history
  or explicit "retired/category-error" annotations.
