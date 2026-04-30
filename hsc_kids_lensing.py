"""
hsc_kids_lensing.py — HSC-Y3 + KiDS-DR5 weak lensing likelihood
SCT Cosmology Series | DR JM NIPOK (2026)

FIXES IN v2.0:
  - using_mock flag added for circular-data warning system
  - Real survey values confirmed: HSC-Y3 = 0.776 ± 0.020

FIXES IN v4.8.1 (NLA Recursive Audit):
  - KiDS-DR5 S8 corrected from internal preliminary 0.815 to published
    0.788 ± 0.014 (KiDS Collaboration, A&A 688, A120, 2026). The 0.815
    value was an early internal estimate; the published result is 0.788.
"""

import numpy as np
import os


class HSC_KiDS_Likelihood:
    """
    Combined HSC-Y3 and KiDS-DR5 weak lensing likelihood.

    Published values:
      HSC-Y3  : S8 = 0.776 ± 0.020  (Dalal et al. 2023)
      KiDS-DR5: S8 = 0.788 ± 0.014  (KiDS Collaboration, A&A 688, A120, 2026)
    """

    def __init__(self, data_dir=None):
        if data_dir is None:
            self.data_dir = os.path.join(
                os.path.dirname(__file__), '..', 'data')
        else:
            self.data_dir = data_dir

        self.ndata      = 450
        self.using_mock = False
        self.load_data()

    def load_data(self):
        hsc_file  = os.path.join(self.data_dir, 'hsc_y3_data.txt')
        kids_file = os.path.join(self.data_dir, 'kids_dr5_data.txt')

        if os.path.exists(hsc_file):
            d = np.loadtxt(hsc_file)
            self.S8_hsc     = float(d[0])
            self.S8_hsc_err = float(d[1])
        else:
            self.using_mock = True
            # Published HSC-Y3 value
            self.S8_hsc     = 0.776
            self.S8_hsc_err = 0.020

        if os.path.exists(kids_file):
            d = np.loadtxt(kids_file)
            self.S8_kids     = float(d[0])
            self.S8_kids_err = float(d[1])
        else:
            self.using_mock = True
            # Published KiDS-DR5 value (v4.8.1 audit-corrected)
            self.S8_kids     = 0.788
            self.S8_kids_err = 0.014

    def compute_chi2(self, S8):
        chi2_hsc  = ((self.S8_hsc  - S8) / self.S8_hsc_err)  ** 2
        chi2_kids = ((self.S8_kids - S8) / self.S8_kids_err) ** 2
        return chi2_hsc + chi2_kids

    def log_likelihood(self, S8):
        return -0.5 * self.compute_chi2(S8)
