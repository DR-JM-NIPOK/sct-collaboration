# sct-collaboration v4.8.3 — change record

## Carried forward from v4.8.2 (A* re-cascade)
A* 5.970→6.173, N_coh 13.51→14.06, f_b unified to 0.162 (cascade-derived); rationale rewritten; CAMB/CLASS/cosmology sector unchanged (R_b=0.2545, c_s²=0.41817, N_eff=2.514).

## New in v4.8.3
### 1. predictions.csv — 32 → 115 (full timeline catalog)
Replaced with the complete 115-prediction catalog from Paper 18 (Confirming Falsifiability), timeline-ordered (Part I 1–25 confirmed/accommodation; Part II 26–115 pending by test era). Current status:
CONFIRMED 2 (P024 R_b, P025 S8) · CONSISTENT 18 · PENDING 90 · CONSISTENCY-REQUIREMENT 2 (P003, P033) · OPEN-DERIVATION 3 (P113–P115). Paper attributions use current numbering.

### 2. Paper cross-references remapped (legacy file-number → current Series/Paper)
Authoritative register map applied to all code/docs (predictions.csv already current). Full map (old file-number → current citation):
| old | new | | old | new |
|---|---|---|---|---|
| 1 | Paper 1 | | 12 | Paper 7 |
| 2 | Paper 2 | | 13 | Paper 12 |
| 3 | Paper 4 | | 14 | Paper 9 |
| 4 | Paper 3 | | 15 | Paper 10 |
| 5 | Paper 5 | | 16 | **Paper 15** (Codified Acoustics) |
| 6 | Paper 11 | | 17 | **Series 2 Paper 1** (Coalescent Parsimony) |
| 7 | Paper 14 | | 18 | Series 2 Paper 2 |
| 8 | Paper 17 | | 19 | Paper 6 |
| 9 | Paper 16 | | 20 | Paper 13 |
| 10 | Paper 18 | | 21–22 | Series 2 Paper 3–4 |

DOI-verified anchors: old "Paper 16" (DOI …10321.29288) = Paper 15; old "Paper 17" (DOI …14355.03366) = Series 2 Paper 1. Ranges like "Papers 1–46" left as series-span references.

## Verified
- predictions.csv = 115 rows; DOI anchors correct; no residual 5.970/0.1675/13.51 in code; corrected RNLA canonical_sweep shows no stale-value or f_b-invariant flags.

## Predictions re-cascade audit (RESOLVED)
- All 115 predictions checked against the A* re-cascade: **none quotes an old A* numeric** (5.970/4.970/13.51/0.1675). Every A*-related prediction is stated qualitatively or as a scaling law (A^(-7/3), Phi_eff = -G*A(r)*M/r, sigma_v/v_bulk, "where A approaches 1") or uses R_b-branch values (R_b, S8, b_IA, N_eff) that the re-cascade does not touch. **No prediction value or status needs changing.**
- The 3 "Paper 17" citations (P086, P094, P097) are **CONFIRMED CORRECT** = current S1 Paper 17 (Constructive Relativity, DOI 10.13140/RG.2.2.23479.79528) - distinct from Series 2 Paper 1 (Coalescent Parsimony, DOI ...14355). No change.
- One A*-derived falsification threshold (P057: M_eff(>50 kpc) > 5e11 M_sun) is a round order-of-magnitude bound, robust to the 3.4%% A* shift - left as is.

## FLAGGED for your review
- **Remap assumption:** the whole site used the legacy file-number scheme (DOI-verified for the two highest-frequency cases, 16 & 17). Spot-check a few of your own if any paper had already been partly updated.
- **Carried-over A* review items** (A_lens P22-equiv, μ_SCT, rotation-fit, tension T023) from v4.8.2 still pending recompute.
- **r_d discrepancy** (separate issue): Paper 18 ledger says r_d=146.8 (provisional) while the site code says 161.4 — not part of this change; reconcile separately.
