/* ============================================================================
 * perturbations_car.c - CAR (late-time) sound speed for CLASS   (v4.8.4)
 * ============================================================================
 *   SCT Cosmology Series | DR JM NIPOK (2026) | License: GPL-3.0
 *
 * RNLA v2.3 CORRECTION (June 2026):
 *   Earlier versions of this file told the user to REPLACE CLASS's
 *   recombination sound speed with cs2 = (1 + R)/3 using the standard density
 *   ratio R = 3*rho_b/(4*rho_g) (~673 at z=0). That is wrong twice over:
 *     (a) it applies the CAR enhancement at recombination - a CATEGORY ERROR
 *         (the CAR enhancement is a LATE-TIME coherent effect that sets S8 and
 *         b_IA; it is NOT the recombination acoustic speed); and
 *     (b) using the density ratio R (~673) in (1 + R)/3 gives cs2 ~ 224 c^2,
 *         i.e. a superluminal speed.
 *
 *   DO NOT modify CLASS's recombination sound speed. The recombination/drag
 *   sound horizon is STANDARD: rs(z_drag) ~ 146.8 Mpc, r_*(z*) ~ 144.4 Mpc,
 *   consistent with DESI-DR2 BAO (147 +/- 1 Mpc).
 *
 *   The CAR coherent sound speed below uses the SCT-DERIVED coherence ratio
 *   R_b(z) = R_B_DERIVED/(1+z) with R_B_DERIVED = 0.2545 (NOT the density
 *   ratio). It belongs in the LATE-TIME growth / weak-lensing / IA sector.
 * ============================================================================ */

/* Recombination sound speed (UNCHANGED - standard CLASS): */
/*   R = (3. * pba->rho_b) / (4. * pba->rho_g);   // density ratio (~673 at z=0) */
/*   pba->cs2 = 1. / (3. * (1. + R));             // standard; -> rs ~ 146.8 Mpc */

/* CAR LATE-TIME coherent sound speed (for the growth/IA sector ONLY): */
/*   R_b = R_B_DERIVED / (1. + z);   // SCT-derived coherence ratio, 0.2545/(1+z) */
/*   cs2_coherent = (1. + R_b) / 3.; // 0.4182 at z=0; -> 1/3 at high z          */

#define R_B_DERIVED_CAR 0.2545   /* Series 2 Paper 1 Section 11.6 (derived) */

/* ============================================================================
 * CAR late-time coherent sound speed (units of c^2).
 * NOTE: this is a LATE-TIME quantity. Do NOT use it for the recombination
 * sound speed or the sound horizon - those stay standard (rs ~ 146.8 Mpc).
 * ============================================================================ */
double CAR_late_time_cs2(double z) {
    double R_b = R_B_DERIVED_CAR / (1.0 + z);  /* SCT-derived, NOT density ratio */
    return (1.0 + R_b) / 3.0;                  /* 0.4182 at z=0; -> 1/3 high-z   */
}

/* ============================================================================
 * Integration note:
 *   - Leave CLASS's recombination sound speed and rs computation UNCHANGED.
 *   - Use CAR_late_time_cs2(z) only in the post-recombination coherent sector
 *     (matter power / S8 normalisation, intrinsic-alignment bias b_IA).
 * ============================================================================ */
