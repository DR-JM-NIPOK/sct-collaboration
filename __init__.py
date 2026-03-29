"""
likelihoods/__init__.py
========================
SCT Cosmology likelihood modules for the Codified Acoustic Relation (CAR).

Available likelihoods
---------------------
    DESIDR2BAOLikelihood        : DESI-DR2 BAO (32 measurements)
    DESY6ThreeTwoPointLikelihood: DES-Y6 3×2pt weak lensing (1800 elements)
    HSCY3WeakLensingLikelihood  : HSC-Y3 cosmic shear (576 elements)
    KiDSDR5WeakLensingLikelihood: KiDS-DR5 weak lensing (912 elements)
    PlanckPR4CMBLikelihood      : Planck PR4 compressed CMB (86 params)

Usage
-----
    from likelihoods import DESIDR2BAOLikelihood
    from sct_core import CAR_predictions

    params = CAR_predictions()
    lik    = DESIDR2BAOLikelihood()
    lnL    = lik.log_like(params)
    chi2   = lik.chi2(params)
"""

from .car_likelihood_base        import CARLikelihoodBase
from .desi_dr2_bao               import DESIDR2BAOLikelihood
from .des_y6_3x2pt               import DESY6ThreeTwoPointLikelihood
from .hsc_y3_wl                  import HSCY3WeakLensingLikelihood
from .kids_dr5_wl                import KiDSDR5WeakLensingLikelihood
from .planck_pr4_cmb             import PlanckPR4CMBLikelihood

__all__ = [
    'CARLikelihoodBase',
    'DESIDR2BAOLikelihood',
    'DESY6ThreeTwoPointLikelihood',
    'HSCY3WeakLensingLikelihood',
    'KiDSDR5WeakLensingLikelihood',
    'PlanckPR4CMBLikelihood',
]
