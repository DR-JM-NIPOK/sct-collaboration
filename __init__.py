"""
likelihoods/__init__.py — Likelihood module initialization
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - Import corrected to match actual filename plank_cmb.py
    (the file has a historical typo — missing 'n' — preserved for
     backwards compatibility; imports updated to match)

FIXES IN v4.8.1 (audit):
  - Re-export DESIDR2BAOLikelihood, DESY6ThreeTwoPointLikelihood, etc.
    so test_likelihoods.py can import them directly from likelihoods.
"""

from .desi_bao import DESI_BAO_Likelihood
from .des_y6_lensing import DES_Y6_Likelihood
from .hsc_kids_lensing import HSC_KiDS_Likelihood
from .plank_cmb import Planck_Likelihood          # file is plank_cmb.py (no 'n')
from .combined_likelihood import CombinedLikelihood, car_loglike

# v4.8.1 audit re-exports
try:
    from .desi_dr2_bao import DESIDR2BAOLikelihood
except Exception:
    DESIDR2BAOLikelihood = None
try:
    from .des_y6_3x2pt import DESY6ThreeTwoPointLikelihood
except Exception:
    DESY6ThreeTwoPointLikelihood = None
try:
    from .hsc_y3_wl import HSCY3WeakLensingLikelihood
except Exception:
    HSCY3WeakLensingLikelihood = None
try:
    from .kids_dr5_wl import KiDSDR5WeakLensingLikelihood
except Exception:
    KiDSDR5WeakLensingLikelihood = None
try:
    from .planck_pr4_cmb import PlanckPR4CMBLikelihood
except Exception:
    PlanckPR4CMBLikelihood = None

__all__ = [
    'DESI_BAO_Likelihood',
    'DES_Y6_Likelihood',
    'HSC_KiDS_Likelihood',
    'Planck_Likelihood',
    'CombinedLikelihood',
    'car_loglike',
]
