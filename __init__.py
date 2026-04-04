"""
likelihoods/__init__.py — Likelihood module initialization
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - Import corrected to match actual filename plank_cmb.py
    (the file has a historical typo — missing 'n' — preserved for
     backwards compatibility; imports updated to match)
"""

from .desi_bao import DESI_BAO_Likelihood
from .des_y6_lensing import DES_Y6_Likelihood
from .hsc_kids_lensing import HSC_KiDS_Likelihood
from .plank_cmb import Planck_Likelihood          # file is plank_cmb.py (no 'n')
from .combined_likelihood import CombinedLikelihood, car_loglike

__all__ = [
    'DESI_BAO_Likelihood',
    'DES_Y6_Likelihood',
    'HSC_KiDS_Likelihood',
    'Planck_Likelihood',
    'CombinedLikelihood',
    'car_loglike',
]
