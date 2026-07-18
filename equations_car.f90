! ============================================================================
! equations_car.f90 — CAR (late-time) sound-speed module for CAMB  (v4.8.4)
! ============================================================================
!
!   SCT Cosmology Series Paper #16 / Paper #17 | DR JM NIPOK (2026)
!   Series 2 Paper 1 DOI: 10.13140/RG.2.2.14355.03366
!   License: GPL-3.0
!
! ----------------------------------------------------------------------------
! v4.8.4 RNLA v2.3 Recursive Audit Correction (June 2026)
! ----------------------------------------------------------------------------
!
! IMPORTANT SCOPE CORRECTION. Earlier versions of this file (v4.8.1) told the
! user to REPLACE CAMB's recombination-epoch sound speed cs2 with the CAR
! formula (1 + R_b)/3 and claimed this yields r_drag ≈ 161.4 Mpc. THAT IS A
! CATEGORY ERROR and is retired:
!
!   • The CAR enhanced sound speed, cs²=(1+R_b)/3 with R_b0=0.2545, is a
!     LATE-TIME coherent-sector quantity. It governs the z≈0 coherence that
!     sets S8 and the intrinsic-alignment bias b_IA. It is NOT the
!     recombination acoustic speed.
!   • At recombination R_b → 0, so (1+R_b)/3 → 1/3 — the pure-radiation speed,
!     which OMITS baryon loading and inflates the horizon. Feeding it into the
!     recombination sound-horizon integral produced the superseded 161.4 Mpc.
!   • The recombination/drag sound horizon is STANDARD, fixed by ω_b and ω_m:
!         r_drag ≈ 146.8 Mpc,   r_*(z*) ≈ 144.4 Mpc   (Planck-consistent).
!
! THEREFORE: do NOT patch CAMB's recombination cs2. Leave CAMB's standard
! photon-baryon acoustic physics untouched; r_drag stays ≈ 146.8 Mpc, which is
! consistent with DESI-DR2 BAO (147 ± 1 Mpc). This module instead provides the
! CAR LATE-TIME coherent sound speed used by the growth / weak-lensing / IA
! sector (S8, b_IA), where the SCT modification actually lives.
!
! ============================================================================
! THEORY
! ============================================================================
!
! Recombination (UNCHANGED — standard CAMB):
!     cs²(z) = 1 / [3 (1 + R(z))],   R(z) = 3 ρ_b / (4 ρ_γ)
!     → r_drag ≈ 146.8 Mpc, r_*(z*) ≈ 144.4 Mpc
!
! Late-time CAR coherent sector (this module):
!     cs²(z) = (1 + R_b(z)) / 3,     R_b(z) = R_B_DERIVED / (1 + z)
!     At z = 0:  cs² = 0.4182  (26% coherent enhancement; sets S8, b_IA)
!     At high z: cs² → 1/3      (coherent enhancement switches off)
!
! R_B_DERIVED = 0.2545 ± 0.032 is fixed by (Series 2 Paper 1 §11.6):
!   (1) SO(3) angular momentum structure of the collision cascade (N=3)
!   (2) QCD phase transition Israel-Darmois junction (Paper 9; 13.6% loss)
! It is a DERIVED constant of the SCT framework, NOT the density ratio R.
!
! ============================================================================
! INTEGRATION INTO CAMB  (late-time / growth sector only)
! ============================================================================
!
! Step 1. Do NOT modify CAMB's recombination cs2 or the tight-coupling
!         sound-speed line. The early-universe acoustic horizon must remain
!         the standard r_drag ≈ 146.8 Mpc.
!
! Step 2. Use CAR_cs2_from_a / CAR_cs2 in the LATE-TIME sector where SCT
!         predicts the coherent enhancement (matter power / S8 normalisation,
!         intrinsic-alignment bias b_IA). These are post-recombination,
!         low-z quantities.
!
! Step 3. Verify the early-time horizon is unchanged:
!           rdrag ≈ 146.8 Mpc   (standard; equals DESI-DR2 within ~0.2σ)
!         and the late-time coherent speed:
!           CAR_cs2(0.0) = 0.4182
!         with: python equations_car_test.py
!
! ============================================================================

module CAR_corrections
  use Precision
  implicit none

  ! Series 2 Paper 1 §11.6 derived constant — DO NOT change without re-deriving
  real(dl), parameter :: R_B_DERIVED = 0.2545_dl

  ! Uncertainty for error propagation in MCMC analyses
  real(dl), parameter :: R_B_UNCERTAINTY = 0.032_dl

contains

  !---------------------------------------------------------------------
  ! CAR_R_b_of_z
  !
  ! SCT-derived baryon-photon coherence ratio at redshift z (LATE-TIME
  ! coherent sector). NOT the standard density ratio R = 3·grho_b/(4·grho_g).
  !---------------------------------------------------------------------
  pure function CAR_R_b_of_z(z) result(Rb)
    real(dl), intent(in) :: z
    real(dl)             :: Rb
    Rb = R_B_DERIVED / (1._dl + z)
  end function CAR_R_b_of_z

  !---------------------------------------------------------------------
  ! CAR_cs2
  !
  ! CAR LATE-TIME coherent sound speed squared (units of c²) at redshift z.
  ! At z = 0:  returns 0.4182  (26% coherent enhancement; sets S8, b_IA)
  ! At high z: returns → 1/3    (enhancement switches off)
  ! WARNING: this is a late-time quantity. Do NOT use it for the
  ! recombination sound horizon (that uses the standard baryon-loaded speed).
  !---------------------------------------------------------------------
  pure function CAR_cs2(z) result(cs2)
    real(dl), intent(in) :: z
    real(dl)             :: cs2
    cs2 = (1._dl + CAR_R_b_of_z(z)) / 3._dl
  end function CAR_cs2

  !---------------------------------------------------------------------
  ! CAR_cs2_from_a
  !
  ! Convenience: CAR late-time coherent sound speed squared vs scale factor.
  !---------------------------------------------------------------------
  pure function CAR_cs2_from_a(a) result(cs2)
    real(dl), intent(in) :: a
    real(dl)             :: cs2, Rb
    Rb  = R_B_DERIVED * a   ! = R_B_DERIVED/(1+z) since a = 1/(1+z)
    cs2 = (1._dl + Rb) / 3._dl
  end function CAR_cs2_from_a

end module CAR_corrections

! ============================================================================
! Example use (LATE-TIME sector only):
!
!     use CAR_corrections, only : CAR_cs2_from_a, R_B_DERIVED
!     ...
!     ! Late-time coherent sound speed for the growth / IA sector:
!     cs2_coherent = CAR_cs2_from_a(a)
!
! NOTE: the recombination sound speed and r_drag are LEFT STANDARD. The CAR
! enhancement is a post-recombination, late-time effect; applying it at
! recombination is the retired 161.4 Mpc category error. Leave CAMB's
! R = 3*grho_b/(4*grho_g) and recombination cs2 unchanged.
! ============================================================================
