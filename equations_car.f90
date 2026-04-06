! ============================================================================
! equations_car.f90 - CAR modification for CAMB
! ============================================================================
! This file contains the modified sound speed calculation for CAMB.
! To use: In your CAMB installation, replace the standard sound speed
! calculation in equations.f90 with the code below.
! ============================================================================

! Standard ΛCDM sound speed (original):
!   R = (3._dl * grho_b) / (4._dl * grho_g)
!   cs2 = 1._dl / (3._dl * (1._dl + R))

! CAR sound speed (modified):
!   R = (3._dl * grho_b) / (4._dl * grho_g)
!   cs2 = (1._dl + R) / 3._dl

! ============================================================================
! Complete modified subroutine for CAMB
! ============================================================================

subroutine CAR_sound_speed(grho_b, grho_g, cs2, R)
  ! NOTE: CAMB uses double precision real(dl) throughout.
  ! Replace 'real' with 'real(dl)' and add 'use Precision' when
  ! integrating into equations.f90. This stub uses real(dl) for correctness.
  use Precision
  implicit none
  real(dl), intent(in) :: grho_b, grho_g
  real(dl), intent(out) :: cs2, R
  ! Baryon-to-photon ratio (same as ΛCDM — computed from density fields, not hardcoded)
  R = (3._dl * grho_b) / (4._dl * grho_g)
  ! CAR sound speed: cs2 = (1+R)/3 replaces standard 1/(3*(1+R))
  cs2 = (1._dl + R) / 3._dl
end subroutine CAR_sound_speed

! ============================================================================
! How to integrate into CAMB:
! 1. Open equations.f90 in your CAMB source directory
! 2. Find the section that calculates the sound speed (search for "cs2")
! 3. Replace the existing calculation with:
!    R = (3._dl * grho_b) / (4._dl * grho_g)
!    cs2 = (1._dl + R) / 3._dl
! 4. Save and recompile CAMB
! ============================================================================