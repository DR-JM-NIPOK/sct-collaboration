! ============================================================================
! equations_car.f90 — CAR modification for CAMB Boltzmann solver  (v4.8.1)
! ============================================================================
!
!   SCT Cosmology Series Paper #16 / Paper #17 | DR JM NIPOK (2026)
!   Paper 17 DOI: 10.13140/RG.2.2.14355.03366
!   License: GPL-3.0
!
! ----------------------------------------------------------------------------
! v4.8.1 NLA Recursive Audit Correction (April 2026)
! ----------------------------------------------------------------------------
!
! THIS FILE WAS CORRECTED in v4.8.1. The previous version of equations_car.f90
! used the standard density-ratio formulation:
!
!     R   = (3 * grho_b) / (4 * grho_g)        ! standard, R ≈ 673 at z=0
!     cs2 = (1._dl + R) / 3._dl                ! gives cs² ≈ 224 at z=0
!
! That formulation produced physically nonsensical results (cs > c at z=0,
! cs² ≈ 224 in units of c²) and was inconsistent with sct_core.py, which
! uses the canonical CAR prescription:
!
!     R_b(z) = R_B_DERIVED / (1 + z)           ! R_B_DERIVED = 0.2545
!     cs²(z) = (1 + R_b(z)) / 3                ! cs² ∈ [1/3, 0.4182]
!
! The correct prescription (matching sct_core.py and Paper 17 v4.8 §11.6)
! is implemented below. After applying this patch and rebuilding CAMB, the
! sound horizon r_drag should be approximately 161.4 Mpc at H0=70.4 km/s/Mpc
! (verified by sct_core.py's compute_r_d_integral).
!
! ============================================================================
! THEORY
! ============================================================================
!
! The Codified Acoustic Relation (CAR):
!
!     Standard ΛCDM:  cs²(z) = 1 / [3 (1 + R(z))]   with R(z) = 3·ρ_b / (4·ρ_γ)
!     CAR ansatz   :  cs²(z) = (1 + R_b(z)) / 3      with R_b(z) = R_B_DERIVED/(1+z)
!
! Critical: R_b ≠ R. The CAR baryon-photon coherence ratio R_b is a DERIVED
! constant of the SCT framework (Paper 17 v4.8 §11.6), NOT the standard
! cosmological density ratio. R_B_DERIVED = 0.2545 ± 0.032 is fixed by:
!   (1) SO(3) angular momentum structure of the collision cascade (N=3)
!   (2) QCD phase transition Israel-Darmois junction (Paper 14, 13.6% loss)
!
! At z_drag ≈ 1060: R_b ≈ 0.0002 → cs² ≈ 1/3 (photon limit)
! At z = 0:         R_b = 0.2545 → cs² = 0.4182 (26% above ΛCDM)
!
! ============================================================================
! INTEGRATION INTO CAMB
! ============================================================================
!
! Step 1. Locate the sound speed calculation in CAMB's equations.f90.
!         (In CAMB 1.5+ it's in subroutine derivs around the line that
!         computes "cs2", typically alongside R = 3*rhob/(4*rhog).)
!
! Step 2. REMOVE the standard ΛCDM lines:
!           R   = 3*rhob/(4*rhog)
!           cs2 = 1._dl/(3._dl*(1._dl + R))
!
! Step 3. REPLACE with the CAR canonical block:
!
!           ! ── CAR canonical sound speed (v4.8.1) ──────────────────
!           ! R_b is the SCT-derived coherence ratio, NOT the density ratio.
!           ! z = 1/a - 1, where a is the scale factor used by CAMB
!           R_B_DERIVED_CAR = 0.2545_dl
!           a_now = 1._dl / (1._dl + z_current)        ! or pass scale factor directly
!           Rb_z  = R_B_DERIVED_CAR * a_now             ! = R_B_DERIVED/(1+z)
!           cs2   = (1._dl + Rb_z) / 3._dl
!           ! ────────────────────────────────────────────────────────
!
! Step 4. Rebuild CAMB with `make`. The patched run should produce:
!           rdrag(H0=70.4) ≈ 161.4 Mpc   [vs ΛCDM 144 Mpc, vs Planck obs 150 Mpc]
!
! Step 5. Verify with: python equations_car_test.py
!
! ============================================================================

module CAR_corrections
  use Precision
  implicit none

  ! Paper 17 v4.8 §11.6 derived constant — DO NOT change without re-deriving
  real(dl), parameter :: R_B_DERIVED = 0.2545_dl

  ! Uncertainty for error propagation in MCMC analyses
  real(dl), parameter :: R_B_UNCERTAINTY = 0.032_dl

contains

  !---------------------------------------------------------------------
  ! CAR_R_b_of_z
  !
  ! Compute the SCT-derived baryon-photon coherence ratio at redshift z.
  ! NOT the standard density ratio (which CAMB calls R = 3·grho_b/(4·grho_g)).
  !---------------------------------------------------------------------
  pure function CAR_R_b_of_z(z) result(Rb)
    real(dl), intent(in) :: z
    real(dl)             :: Rb
    Rb = R_B_DERIVED / (1._dl + z)
  end function CAR_R_b_of_z

  !---------------------------------------------------------------------
  ! CAR_cs2
  !
  ! CAR sound speed squared (in units of c²) at redshift z.
  ! At z_drag ≈ 1060: returns ≈ 1/3 + 8e-5 (essentially photon limit)
  ! At z = 0:         returns 0.4182  (26% above ΛCDM photon limit)
  !---------------------------------------------------------------------
  pure function CAR_cs2(z) result(cs2)
    real(dl), intent(in) :: z
    real(dl)             :: cs2
    cs2 = (1._dl + CAR_R_b_of_z(z)) / 3._dl
  end function CAR_cs2

  !---------------------------------------------------------------------
  ! CAR_cs2_from_a
  !
  ! Convenience: CAR sound speed squared as a function of scale factor.
  ! For CAMB integration, scale factor is often more accessible than z.
  !---------------------------------------------------------------------
  pure function CAR_cs2_from_a(a) result(cs2)
    real(dl), intent(in) :: a
    real(dl)             :: cs2, Rb
    Rb  = R_B_DERIVED * a   ! = R_B_DERIVED/(1+z) since a = 1/(1+z)
    cs2 = (1._dl + Rb) / 3._dl
  end function CAR_cs2_from_a

end module CAR_corrections

! ============================================================================
! Example replacement block for CAMB equations.f90:
!
!     use CAR_corrections, only : CAR_cs2_from_a, R_B_DERIVED
!     ...
!     ! ── CAR canonical (replaces ΛCDM cs2 line) ──────────────────────────
!     cs2 = CAR_cs2_from_a(a)
!     ! ────────────────────────────────────────────────────────────────────
!
! NOTE: The standard density-ratio "R" symbol used internally by CAMB
! (R = 3*grho_b/(4*grho_g)) is NO LONGER needed for cs2 in the CAR
! prescription. It may still be used elsewhere in CAMB (e.g., for
! photon-baryon momentum coupling); leave those uses unchanged. Only
! the cs2 calculation should switch to the CAR formula.
! ============================================================================
