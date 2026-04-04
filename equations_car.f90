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
  implicit none
  real, intent(in) :: grho_b, grho_g
  real, intent(out) :: cs2, R
  ! Baryon-to-photon ratio (same as ΛCDM)
  R = (3.0 * grho_b) / (4.0 * grho_g)
  ! CAR sound speed (modified)
  cs2 = (1.0 + R) / 3.0
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