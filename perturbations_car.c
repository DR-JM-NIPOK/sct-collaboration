/* ============================================================================
 * perturbations_car.c - CAR modification for CLASS
 * ============================================================================
 * This file contains the modified sound speed calculation for CLASS.
 * To use: In your CLASS installation, replace the standard sound speed
 * calculation in perturbations.c with the code below.
 * ============================================================================ */

/* Standard ΛCDM sound speed (original): */
/*   R = (3. * pba->rho_b) / (4. * pba->rho_g); */
/*   pba->cs2 = 1. / (3. * (1. + R)); */

/* CAR sound speed (modified): */
/*   R = (3. * pba->rho_b) / (4. * pba->rho_g); */
/*   pba->cs2 = (1. + R) / 3.; */

/* ============================================================================
 * Complete modified function for CLASS
 * ============================================================================ */

void CAR_sound_speed(struct background *pba, double *cs2, double *R) {
    /* Baryon-to-photon ratio (same as ΛCDM) */
    *R = (3.0 * pba->rho_b) / (4.0 * pba->rho_g);
    /* CAR sound speed (modified) */
    *cs2 = (1.0 + *R) / 3.0;
}

/* ============================================================================
 * How to integrate into CLASS:
 * 1. Open perturbations.c in your CLASS source directory
 * 2. Find the section that calculates the sound speed (search for "cs2" or "pba->cs2")
 * 3. Replace the existing calculation with:
 *    R = (3. * pba->rho_b) / (4. * pba->rho_g);
 *    pba->cs2 = (1. + R) / 3.;
 * 4. Save and recompile CLASS
 * ============================================================================ */