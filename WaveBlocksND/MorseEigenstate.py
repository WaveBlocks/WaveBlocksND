"""The WaveBlocks Project

This file contains the Morse eigenstates.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import arange, zeros, exp, sqrt, floor, linspace, complexfloating, sum, zeros_like, pi, log
from scipy.special import gamma, eval_genlaguerre

from WaveBlocksND.Wavepacket import Wavepacket

__all__ = ["MorseEigenstate"]


class MorseEigenstate(Wavepacket):
    r"""This class implements the abstract :py:class:`Wavepacket` interface
    and contains code for the Morse eigenstates.
    """

    def __init__(self, eps, beta, V0):
        r"""
        A Morse eigenstate.

        :param eps: The semiclassical scaling parameter.
        :param beta: Parameter describing the Morse potential.
        :param V0: Parameter describing the Morse potential.
        """
        self._dimension = 1
        self._number_components = 1

        self._eps = eps

        self._beta = beta
        self._V0 = V0

        self._nu = sqrt(8 * V0 / (beta**2 * eps**4))


    def get_nu(self):
        r"""The value :math:`nu = \sqrt{\frac{8 V_0}{\beta^2 \varepsilon^4}}`.
        """
        return self._nu


    def get_max_levels(self):
        r"""The maximal number of eigenstates possible with the given values of :math:`\beta` and :math:`V_0`.
        """
        return int(floor((self._nu - 1) / 2))


    def _Nn(self, n):
        r"""The normalization constant :math:`N_n`.

        :param n: Number of the eigenstate.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`n` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        return sqrt(self._beta * (self._nu - 2 * n - 1) * gamma(n + 1) / gamma(self._nu - n))


    def _mun(self, n, x):
        r"""Evaluate the :math:`n`-th eigenstate :math:`\mu_n` by direct computation via the analytic formula.

        :param n: Number of the eigenstate.
        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`n` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        sn = 1 / 2 * (self._nu - 2 * n - 1)
        z = self._nu * exp(-self._beta * x.reshape(-1))
        mun = self.Nn(n) * exp(-z / 2) * z**sn * eval_genlaguerre(n, 2 * sn, z)
        return mun


    def _mu0(self, x):
        r"""Evaluate the groundstate :math:`\mu_0` by direct computation
        via an improved (but approximate) analytic formula.

        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.
        """
        # Improved version for small epsilon
        delta = ((self._nu - 1) + 1 / (12 * (self._nu - 1) - 1 / (10 * (self._nu - 1))))
        r2 = sqrt(self._beta * delta / (exp(1) * self._nu))
        r4 = ((self._nu - 1) / (2 * pi)) ** 0.25
        ex = log(delta) - 1 - log(self._nu) + (self._nu - 1) / self._nu * self._beta * x + exp(-self._beta * x)
        return r2 * r4 * exp(-self._nu / 2 * ex)


    def evaluate_direct(self, K, x):
        r"""Evaluate the first :math:`K` eigenstates :math:`\mu_n` by direct computation via the analytic formula.

        :param K: List of the eigenstates to compute.
        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`n` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        Kmax = K.shape[0]
        g = x.shape[1]
        B = zeros((Kmax, g), dtype=complexfloating)

        for n in K:
            B[n, :] = self.mu(n, x)
        return B


    def evaluate_recursive(self, K, x):
        r"""Evaluate the first :math:`K` eigenstates :math:`\mu_n` by a recursive scheme, starting
        with an approximation to the groundstate :math:`\mu_0`.

        :param K: List of the eigenstates to compute.
        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.
        """
        Kmax = K.shape[0]
        g = x.shape[1]
        B = zeros((Kmax, g), dtype=complexfloating)

        # Groundstate
        B[0, :] = self._mu0(x)

        # Recursion
        for n in range(0, Kmax - 1):
            rn = sqrt((n + 1) * (self._nu - n - 1))
            ln = sqrt(n * (self._nu - n))
            sn = 1 / 2 * (self._nu - 2 * n - 1)
            pf1 =  1 / rn * sqrt((sn - 1) / sn      ) * (2 * sn)     / (2 * sn + 1)
            pf2 = ln / rn * sqrt((sn - 1) / (sn + 1)) * (2 * sn - 1) / (2 * sn + 1)
            pf1 *= ((4 * sn**2 - 1) * exp(self._beta * x) / self._nu - self._nu)
            # Note: Index wraps around for n = 0
            B[n + 1, :] = pf1 * B[n, :] - pf2 * B[n - 1, :]

        return B


    def En(self, n):
        r"""Compute the energy of the :math:`n`-th eigenstate :math:`\mu_n` by the analytic formula.

        :param n: Number of the eigenstate.
        """
        return -1 / 2 * (sqrt(2 * self._V0) - self._eps**2 * abs(self._beta) * (n + 1 / 2))**2
