# sct-collaboration v4.8.5 — change record (register update)

Source: in-session verification note `SCT_VERIFICATION_NOTE_CAMB_SPARC_20260710.md`
(CAMB 1.6.6, pip/Linux sandbox; Planck 2018 binned TT spectrum, 83 band powers;
SPARC Lelli+2016c tables, Zenodo 16284118 mirror). This is a **documentation and
register-consistency release** — no theoretical or numeric constant changes.
Base: v4.8.4.

## What changed

### N_eff = 2.514 quantified against Planck (new tensions.csv T031)
The verification note ran a first-order, stock-mapping χ² test: substituting
N_eff = 2.514 into standard ΛCDM physics at fixed θ* against the Planck 2018
binned TT spectrum gives **Δχ² ≈ +658** over 83 band powers (H0 dragged to
≈61.3 km/s/Mpc). This decisively excludes the *naive* substitution and
quantifies, for the first time in-repository, the 2.8σ tension already noted
qualitatively in `predictions.csv` P088 (SCT 2.514 vs Planck 2.99 ± 0.17).
Logged as `tensions.csv` **T031** (status: Open/PARTIAL) and folded into
P088's notes.

### Register-consistency flag — which era does N_eff = 2.514 describe?
The note surfaces a real internal tension in how the register currently reads:
the H0 ≈ 66–67 note (global, from θ* + standard sound horizon) only reproduces
under **standard** N_eff = 3.044 (θ*+standard EOS → H0 = 65.93). Combining that
same EOS with **N_eff = 2.514** instead gives H0 = 60.15 — inconsistent with the
register's own H0 claim. Under Plasma Equivalence (standard recombination-era
physics), the natural reading is that N_eff = 2.514 is a derived *late-time*
quantity, not a recombination-era one — but the existing register language
frames 2.514 as literally what CMB-S4 will measure. Both readings cannot hold
simultaneously; `predictions.csv` P088 and `README.md` (root and
`paper17-v4-deriver-Rb/`) now flag this explicitly and require the CAR-patched
Boltzmann run to resolve which observable N_eff = 2.514 actually is.

### Audit-trail cross-reference (informational only, no value changed)
r_drag evaluated at N_eff = 2.514 (149.86 Mpc) numerically lands close to the
retired 149.1 Mpc r_d figure from pre-v4.8.4 drafts. The note flags this as a
coincidence worth remembering during future audits — it is **not** a revival of
the retracted value. r_d remains 146.8 Mpc (v4.8.4, RNLA v2.3), unchanged here.

### README.md stale-statement fix
Root `README.md` previously stated "Paper 15 will be rewritten to reflect
this [r_d fix]" in future tense. Cross-checked against the author's current
Paper 15 draft (v2.22, dated March 27 2026 in the manuscript masthead, file
last modified 2026-06-02): the paper already states r_d = 146.8 ± 5 Mpc with
149.1 and 161.4 Mpc explicitly marked SUPERSEDED, using the same reasoning as
the v4.8.4 RNLA v2.3 fix. The README line is corrected to reflect that this is
already done, not pending.

## Files changed
- `tensions.csv` x2 (root, `paper17-v4-deriver-Rb/`): new row T031.
- `predictions.csv` x2: P088 notes expanded; status -> "PENDING (PARTIAL
  sub-check, see notes)".
- `sct_core.py` x10: VERSION HISTORY entry; validation/compact report
  headers and footers relabeled v4.8.5; new "REGISTER UPDATE - v4.8.5" note
  block added alongside (not replacing) the existing v4.8.4 AUDIT NOTE block.
  All printed numeric values verified byte-identical to v4.8.4.
- `README.md` (root): N_eff table row footnoted; stale "will be rewritten"
  line corrected; header bumped to v4.8.5.
- `paper17-v4-deriver-Rb/README.md`: same N_eff register-consistency caveat
  added near its Key Constants table (its own "Version 3.0" banner and
  R_b=0.260-legacy-warning content are untouched — out of scope here).
- `CHANGELOG.md` x10: this entry appended (append-only history), in both the
  longer root variant and the shorter 9-branch variant.

## Out of scope
The verification note's SPARC/RAR full-sample recovery (its item W224) and
CMB tilt-running check (W204) concern papers S1_09 and S1_04. Neither paper is
part of this code repository, and no repository files reference them, so
nothing here was changed for those items.

## Unchanged (verified still canonical)
A* = 6.173, N_coh = 14.06, f_b = 0.162; R_b = 0.2545, c_s² = 0.41817,
b_IA = 1.0848, r_d = 146.8 Mpc, r_*(z*) = 144.4 Mpc, H0 (global) ≈ 66.3,
N_eff = 2.514 (value unchanged — only its status note updated). Predictions
catalog count unchanged at 115.

## Verification
- `python3 sct_core.py --validate` and `python3 sct_core.py` (compact report):
  numeric outputs identical to v4.8.4; only version labels and the new audit
  note block differ.
- `python3 -m py_compile sct_core.py`: passes on all 10 branch copies
  (byte-identical, single md5 across all copies).
- `tensions.csv` / `predictions.csv`: re-parsed with Python's `csv` module
  after editing to confirm well-formed rows; root and `paper17-v4-deriver-Rb/`
  copies confirmed byte-identical to each other post-edit.
