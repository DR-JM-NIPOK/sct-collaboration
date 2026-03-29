"""
likelihoods/car_likelihood_base.py
===================================
Abstract base class for all CAR likelihood modules.

Each survey likelihood inherits from CARLikelihoodBase and implements:
    - load_data()   : load data vector and covariance matrix
    - theory()      : compute theoretical prediction given CAR parameters
    - log_like()    : return log-likelihood scalar

All likelihoods accept the same CAR parameter dictionary produced by
sct_core.CAR_predictions().

Author : DR JM NIPOK
License: GPL-3.0
"""

import abc
import numpy as np
from typing import Optional
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from sct_core import CAR_predictions, BBN_OMEGA_B_H2, PLANCK_OMEGA_M


class CARLikelihoodBase(abc.ABC):
    """Abstract base class for CAR likelihood modules."""

    # Subclasses set these at class level
    name: str = "base"
    n_data: int = 0

    def __init__(self, data_path: Optional[str] = None, verbose: bool = False):
        self.data_path = data_path
        self.verbose   = verbose
        self.data_vec  = None   # shape (n_data,)
        self.cov       = None   # shape (n_data, n_data)
        self.inv_cov   = None   # shape (n_data, n_data)
        self._loaded   = False

    def _ensure_loaded(self):
        if not self._loaded:
            self.load_data()
            if self.cov is not None:
                self.inv_cov = np.linalg.inv(self.cov)
            self._loaded = True

    @abc.abstractmethod
    def load_data(self) -> None:
        """Load data vector and covariance matrix into self.data_vec and self.cov."""

    @abc.abstractmethod
    def theory(self, params: dict) -> np.ndarray:
        """
        Compute theoretical prediction vector for given CAR parameters.

        Parameters
        ----------
        params : dict
            Output of sct_core.CAR_predictions().

        Returns
        -------
        np.ndarray
            Theory vector, same shape as self.data_vec.
        """

    def log_like(self, params: dict) -> float:
        """
        Compute Gaussian log-likelihood:
            ln L = -½ Δᵀ C⁻¹ Δ    where Δ = data − theory

        Parameters
        ----------
        params : dict
            Output of sct_core.CAR_predictions().

        Returns
        -------
        float
            Log-likelihood scalar.
        """
        self._ensure_loaded()
        theory_vec = self.theory(params)
        delta = self.data_vec - theory_vec
        chi2 = float(delta @ self.inv_cov @ delta)
        if self.verbose:
            n = len(delta)
            print(f"  [{self.name}] χ²={chi2:.2f}, dof={n}, χ²/dof={chi2/n:.3f}")
        return -0.5 * chi2

    def chi2(self, params: dict) -> float:
        """Return χ² = -2 ln L."""
        return -2.0 * self.log_like(params)

    def chi2_per_dof(self, params: dict) -> float:
        """Return χ²/dof."""
        self._ensure_loaded()
        return self.chi2(params) / self.n_data

    def __repr__(self):
        return f"{self.__class__.__name__}(name={self.name!r}, n_data={self.n_data})"
