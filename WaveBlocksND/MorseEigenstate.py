"""The WaveBlocks Project

This file contains the Morse eigenstates.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, exp, sqrt, floor, complexfloating, pi, log, array, arange, sum
from scipy.linalg import norm
from scipy.special import gamma, eval_genlaguerre

from WaveBlocksND.Wavepacket import Wavepacket
from WaveBlocksND.HyperCubicShape import HyperCubicShape

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
        self._Kmax = int(floor((self._nu - 1) / 2))


        # Default basis shape
        self._basis_shape = HyperCubicShape([1])

        # A Groundstate
        self._coefficients = zeros((self._basis_shape.get_basis_size(), 1), dtype=complexfloating)

        # Cache basis sizes
        self._basis_size = self._basis_shape.get_basis_size()


    def get_eps(self):
        r"""Retrieve the semi-classical scaling parameter :math:`\varepsilon` of the wavepacket.

        :return: The value of :math:`\varepsilon`.
        """
        return self._eps


    def get_nu(self):
        r"""The value :math:`nu = \sqrt{\frac{8 V_0}{\beta^2 \varepsilon^4}}`.
        """
        return self._nu


    def get_max_levels(self):
        r"""The maximal number of eigenstates possible with the given values of :math:`\beta` and :math:`V_0`.
        """
        return self._Kmax


    def _Nn(self, n):
        r"""The normalization constant :math:`N_n`.

        :param n: Number of the eigenstate.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`n` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        return sqrt(self._beta * (self._nu - 2 * n - 1) * gamma(n + 1) / gamma(self._nu - n))


    def _resize_coefficient_storage(self, bs_old, bs_new):
        r"""
        """
        bso = bs_old.get_basis_size()
        bsn = bs_new.get_basis_size()

        # Find the intersection of K and K'
        # Optimization: iterate over smaller set
        if bso <= bsn:
            insec = [k for k in bs_old if k in bs_new]
        elif bso > bsn:
            insec = [k for k in bs_new if k in bs_old]
        # TODO: Consider making this part of the BasisShape interface
        # TODO: Consider implementing set operations for basis shapes

        # Construct the index mapping
        i = array([bs_old[k] for k in insec])
        j = array([bs_new[k] for k in insec])

        # Copy over the data
        cnew = zeros((bsn, 1), dtype=complexfloating)
        cnew[j] = self._coefficients[i]
        self._coefficients = cnew


    def get_basis_shapes(self):
        r"""Retrieve the basis shape :math:`\mathfrak{K}`.

        :return: The basis shape.
        """
        return self._basis_shape


    def set_basis_shapes(self, basis_shape):
        r"""Set the basis shape :math:`\mathfrak{K}`.

        :param basis_shape: The basis shape.
        :type basis_shape: A subclass of :py:class:`BasisShape`.
        """
        basis_size = basis_shape.get_basis_size()
        if basis_size > self._Kmax:
            raise ValueError("Basis shape too large (size: {}) because only {} eigenstates exist.".format(basis_size, self._Kmax))

        # Adapt the coefficient storage vectors
        self._resize_coefficient_storage(self._basis_shape, basis_shape)
        # Set the new basis shape for the given component
        self._basis_shape = basis_shape
        # And update the caches information
        self._basis_size = basis_size


    def set_coefficient(self, index, value):
        r"""Set a single coefficient :math:`c^i` of :math:`\mu_i`.

        :param index: The multi-index :math:`k` of the coefficient :math:`c^i_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :param value: The new value of the coefficient :math:`c^i_k`.
        :raise: :py:class:`ValueError` For invalid indices :math:`i`.
        """
        if index not in self._basis_shape:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shape[index]
        self._coefficients[key] = value


    def get_coefficient(self, index):
        r"""Retrieve a single coefficient :math:`c^i` of :math:`\mu_i`.

        :param index: The multi-index :math:`k` of the coefficient :math:`c^i_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :return: A single complex number.
        :raise: :py:class:`ValueError` For invalid indices :math:`i`.
        """
        if index not in self._basis_shape:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shape[index]
        return self._coefficients[key]


    def set_coefficients(self, values):
        r"""Update all the coefficients :math:`c` of :math:`\mu`.

        Note: this method copies the data arrays.

        :param values: The new values of the coefficients :math:`c^i` of :math:`\Phi_i`.
        :type values: An ndarray of suitable shape or a list of ndarrays.
        """
        bs = self._basis_size
        self._coefficients = values.copy().reshape((bs, 1))


    def get_coefficients(self):
        r"""Returns the coefficients :math:`c` of :math:`\mu`

        Note: this method copies the data arrays.

        :return: A single ndarray with the coefficients of the given component or
                 a list containing the ndarrays for each component. Each ndarray
                 is two-dimensional with a shape of :math:`(|\mathfrak{K}_i|, 1)`.
        """
        return self._coefficients.copy()


    def get_coefficient_vector(self):
        r"""Retrieve the coefficients :math:`c` of :math:`\mu`.

        .. warning:: This function does *not* copy the input data!
                     This is for efficiency as this routine is used in the innermost loops.

        :return: The coefficients :math:`c` of :math:`\mu`.
        """
        return self._coefficients


    def set_coefficient_vector(self, vector):
        """Set the coefficients of :math:`\mu`.

        .. warning:: This function does *not* copy the input data!
                     This is for efficiency as this routine is used in the innermost loops.

        :param vector: The coefficients of all components as a single long column vector.
        :type vector: A two-dimensional ndarray of appropriate shape.
        """
        self._coefficients = vector


    def _evaluate_mun(self, n, x):
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


    def _evaluate_mu0(self, x):
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


    def _evaluate_basis_at_direct(self, x):
        r"""Evaluate the eigenstates :math:`\mu_n` by direct computation via the analytic formula.

        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`n` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        g = x.shape[1]
        B = zeros((self._basis_size, g), dtype=complexfloating)

        for n in range(self._basis_size):
            B[n, :] = self._evaluate_mun(n, x)
        return B


    def evaluate_basis_at(self, x):
        r"""Evaluate the eigenstates :math:`\mu_n` by a recursive scheme, starting
        with an approximation to the groundstate :math:`\mu_0`.

        :param x: Array of grid nodes to evaluate :math:`\mu_n` on.
        """
        g = x.shape[1]
        B = zeros((self._basis_size, g), dtype=complexfloating)

        # Groundstate
        B[0, :] = self._evaluate_mu0(x)

        # Recursion
        for n in range(0, self._basis_size - 1):
            rn = sqrt((n + 1) * (self._nu - n - 1))
            ln = sqrt(n * (self._nu - n))
            sn = 1 / 2 * (self._nu - 2 * n - 1)
            pf1 =  1 / rn * sqrt((sn - 1) / sn      ) * (2 * sn)     / (2 * sn + 1)
            pf2 = ln / rn * sqrt((sn - 1) / (sn + 1)) * (2 * sn - 1) / (2 * sn + 1)
            pf1 *= ((4 * sn**2 - 1) * exp(self._beta * x) / self._nu - self._nu)
            # Note: Index wraps around for n = 0
            B[n + 1, :] = pf1 * B[n, :] - pf2 * B[n - 1, :]

        return B


    def evaluate_at(self, x):
        r"""Evaluate the Hagedorn wavepacket :math:`\Psi` at the given nodes :math:`\gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
                          (Defaults to ``None`` for evaluating all components.)
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :return: A list of arrays or a single array containing the values of the :math:`\Phi_i` at the nodes :math:`\gamma`.
        """
        B = self.evaluate_basis_at(x)
        mu = sum(self._coefficients * B, axis=0).reshape(1, -1)
        return mu


    def norm(self):
        r"""Calculate the :math:`L^2` norm :math:`\mu\Psi|\mu\rangle` of the wavepacket :math:`\mu`.

        :return: The norm of :math:`\mu`.
        """
        return norm(self._coefficients)


    def energy_levels(self):
        r"""Compute the energy of the :math:`n`-th eigenstate :math:`\mu_n` by the analytic formula.

        :param n: Number of the eigenstate.
        """
        K = arange(self._basis_size)
        return -1 / 2 * (sqrt(2 * self._V0) - self._eps**2 * abs(self._beta) * (K + 1 / 2))**2
