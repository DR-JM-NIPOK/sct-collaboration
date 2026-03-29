"""
chains/run_polychord.py
========================
PolyChord nested sampling runner for the CAR framework.

Implements the full joint likelihood combining all four datasets:
    DESI-DR2 BAO + Planck PR4 CMB + DES-Y6 3×2pt + HSC-Y3 + KiDS-DR5

Parameter space:
    CAR  : Ω_m ∈ [0.1, 0.5],  R_b ∈ [0.24, 0.28]  (2 free parameters)
    ΛCDM : Ω_b h², Ω_c h², 100θ_s, τ, ln(10¹⁰A_s), n_s  (6 cosmological)
           + up to 42 nuisance parameters (IA, photo-z, baryonic feedback)

Author : DR JM NIPOK | License: GPL-3.0
"""

import argparse
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sct_core import CAR_predictions, BBN_OMEGA_B_H2, PLANCK_OMEGA_GAM_H2


def build_combined_likelihood(data: str, verbose: bool = False):
    """
    Construct the joint likelihood for all datasets.

    Parameters
    ----------
    data : str
        Which datasets to include: 'combined', 'desi', 'des', 'hsc', 'kids'
    verbose : bool

    Returns
    -------
    callable
        log_like(theta) function accepting parameter vector theta.
    """
    from likelihoods import (
        DESIDR2BAOLikelihood,
        DESY6ThreeTwoPointLikelihood,
        HSCY3WeakLensingLikelihood,
        KiDSDR5WeakLensingLikelihood,
        PlanckPR4CMBLikelihood,
    )

    likelihoods = []
    if data in ('combined', 'desi'):
        likelihoods.append(DESIDR2BAOLikelihood(verbose=verbose))
    if data in ('combined', 'planck'):
        likelihoods.append(PlanckPR4CMBLikelihood(verbose=verbose))
    if data in ('combined', 'des'):
        likelihoods.append(DESY6ThreeTwoPointLikelihood(verbose=verbose))
    if data in ('combined', 'hsc', 'hsc_kids'):
        likelihoods.append(HSCY3WeakLensingLikelihood(verbose=verbose))
    if data in ('combined', 'kids', 'hsc_kids'):
        likelihoods.append(KiDSDR5WeakLensingLikelihood(verbose=verbose))

    def log_like_car(theta):
        """
        Log-likelihood for CAR model.
        theta = [Omega_m, R_b]
        """
        Omega_m, R_b = theta
        # Reconstruct Omega_b_h2 from R_b
        Omega_b_h2 = R_b * 3.0 * PLANCK_OMEGA_GAM_H2 / 4.0
        try:
            params = CAR_predictions(
                Omega_b_h2 = Omega_b_h2,
                Omega_m    = Omega_m,
            )
            params['Omega_m'] = Omega_m
            params['R_b']     = R_b
        except Exception:
            return -1e30

        total_lnL = sum(lik.log_like(params) for lik in likelihoods)
        return total_lnL

    return log_like_car, len(likelihoods)


def prior_transform_car(cube):
    """
    Map unit hypercube to CAR prior.
    theta[0] = Omega_m  ~ Uniform[0.1, 0.5]
    theta[1] = R_b      ~ Gaussian(0.260, 0.002)  [BBN-anchored]
    """
    from scipy.stats import norm
    theta = np.zeros(2)
    theta[0] = 0.1 + cube[0] * 0.4          # Omega_m: uniform [0.1, 0.5]
    theta[1] = norm.ppf(cube[1], 0.260, 0.002)  # R_b: Gaussian BBN prior
    return theta


def run_polychord_car(model: str, data: str, output_dir: str,
                      n_live: int = 500, n_repeats: int = 10,
                      precision: float = 0.01, n_threads: int = 4):
    """Run PolyChord nested sampling for the CAR or ΛCDM model."""
    try:
        import pypolychord
        from pypolychord.settings import PolyChordSettings
    except ImportError:
        print("ERROR: pypolychord not installed. See SETUP_INSTRUCTIONS.md.")
        print("       Install with: pip install pypolychord  (or build from source)")
        sys.exit(1)

    if model == 'car':
        nDims    = 2
        nDerived = 4   # r_d, H0, S8, IA_bias
        log_like_fn, n_liks = build_combined_likelihood(data)

        def log_like(theta):
            lnL = log_like_fn(theta)
            Omega_m, R_b = theta
            Omega_b_h2 = R_b * 3.0 * PLANCK_OMEGA_GAM_H2 / 4.0
            params = CAR_predictions(Omega_b_h2=Omega_b_h2, Omega_m=Omega_m)
            phi = [params['r_d_Mpc'], params['H0'], params['S8'], params['IA_bias']]
            return lnL, phi

        def prior(cube):
            return prior_transform_car(cube)

        param_names = ['Omega_m', 'R_b', 'r_d', 'H0', 'S8', 'IA_bias']

    else:
        print(f"Model '{model}' requires full ΛCDM CAMB integration.")
        print("See chains/config_lcdm.ini for CosmoMC configuration.")
        return

    settings = PolyChordSettings(nDims, nDerived)
    settings.file_root         = os.path.join(output_dir, f'{model}_{data}')
    settings.n_live            = n_live
    settings.num_repeats       = n_repeats
    settings.precision_criterion = precision
    settings.nprior            = 500
    settings.do_clustering     = True
    settings.read_resume       = False
    settings.write_resume      = True
    settings.feedback          = 1

    print(f"  Running PolyChord: nDims={nDims}, n_live={n_live}, "
          f"n_repeats={n_repeats}")
    pypolychord.run_polychord(log_like, nDims, nDerived, settings, prior,
                               param_names=param_names)
    print(f"  PolyChord complete. Output in: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run PolyChord for CAR model")
    parser.add_argument("--model",     default="car",      choices=["car","lcdm"])
    parser.add_argument("--data",      default="combined")
    parser.add_argument("--output",    default="./output/chains")
    parser.add_argument("--nlive",     type=int,   default=500)
    parser.add_argument("--repeats",   type=int,   default=10)
    parser.add_argument("--precision", type=float, default=0.01)
    parser.add_argument("--threads",   type=int,   default=4)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    run_polychord_car(
        model      = args.model,
        data       = args.data,
        output_dir = args.output,
        n_live     = args.nlive,
        n_repeats  = args.repeats,
        precision  = args.precision,
        n_threads  = args.threads,
    )
