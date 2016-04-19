"""The WaveBlocks Project

This file contains the Morse eigenstates.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, exp, sqrt, floor, complexfloating, pi, log, array, arange, sum, atleast_2d
from scipy.linalg import norm
from scipy.special import gamma, eval_genlaguerre

from WaveBlocksND.Wavepacket import Wavepacket
from WaveBlocksND.AbstractGrid import AbstractGrid
from WaveBlocksND.GridWrapper import GridWrapper
from WaveBlocksND.HyperCubicShape import HyperCubicShape

__all__ = ["MorseEigenstate"]


class MorseEigenstate(Wavepacket):
    r"""This class implements the abstract :py:class:`Wavepacket` interface
    and contains code for the Morse eigenstates.
    """

    def __init__(self, eps, beta, V0):
        r"""
        An eigenstate :math:`\mu` of the Morse potential:

        .. math:: M(x) := V_0 \left(e^{-2 \beta x} - 2 e^{-\beta x}\right)

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
        r"""The value:

        .. math:: \nu = \sqrt{\frac{8 V_0}{\beta^2 \varepsilon^4}}

        :return: The value of :math:`\nu`.
        """
        return self._nu


    def get_max_levels(self):
        r"""The maximal number of eigenstates possible with the given values of :math:`\beta`, :math:`V_0` and :math:`\varepsilon`:

        .. math:: n_{\max} = \left\lfloor \frac{\nu - 1}{2} \right\rfloor
        """
        return self._Kmax


    def _Nn(self, n):
        r"""The normalization constant:

        .. math:: N_n := \sqrt{\beta \frac{(\nu-2n-1) \Gamma(n+1)}{\Gamma(\nu-n)} }

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
        r"""Set a single coefficient :math:`c_k` of :math:`\mu_k`.

        :param index: The multi-index :math:`k` of the coefficient :math:`c_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :param value: The new value of the coefficient :math:`c_k`.
        :raise: :py:class:`ValueError` For invalid indices :math:`k`.
        """
        if index not in self._basis_shape:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shape[index]
        self._coefficients[key] = value


    def get_coefficient(self, index):
        r"""Retrieve a single coefficient :math:`c_k` of :math:`\mu_k`.

        :param index: The multi-index :math:`k` of the coefficient :math:`c_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :return: A single complex number.
        :raise: :py:class:`ValueError` For invalid indices :math:`k`.
        """
        if index not in self._basis_shape:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shape[index]
        return self._coefficients[key]


    def set_coefficients(self, values):
        r"""Update all the coefficients :math:`c` of :math:`\mu`.

        .. note:: This method copies the data arrays.

        :param values: The new values of the coefficients :math:`c` of :math:`\mu`.
        :type values: An ndarray of shape :math:`(|\mathfrak{K}|, 1)`.
        """
        bs = self._basis_size
        self._coefficients = values.copy().reshape((bs, 1))


    def get_coefficients(self):
        r"""Returns the coefficients :math:`c` of :math:`\mu`

        .. note:: This method copies the data arrays.

        :return: A single ndarray with the coefficients. The ndarray
                 is two-dimensional with a shape of :math:`(|\mathfrak{K}|, 1)`.
        """
        return self._coefficients.copy()


    def get_coefficient_vector(self):
        r"""Retrieve the coefficients :math:`c` of :math:`\mu`.

        .. note:: This function does *not* copy the input data!
                  This is for efficiency as this routine is used in the innermost loops.

        :return: The coefficients :math:`c` of :math:`\mu` as an ndarray of shape :math:`(|\mathfrak{K}|, 1)`.
        """
        return self._coefficients


    def set_coefficient_vector(self, vector):
        """Set the coefficients of :math:`\mu`.

        .. note:: This function does *not* copy the input data!
                  This is for efficiency as this routine is used in the innermost loops.

        :param vector: The coefficients as a single long column vector.
        :type vector: A two-dimensional ndarray of shape :math:`(|\mathfrak{K}|, 1)`.
        """
        self._coefficients = vector


    def _evaluate_mun(self, n, x):
        r"""Evaluate the :math:`n`-th eigenstate :math:`\mu_n` by direct computation via the analytic formula:

        .. math:: \mu_n(x) := N_n \, e^{-\frac{z}{2}} \, z^{s_n} \, \mathrm{L}_n^{2{s_n}}(z)

        where:

        .. math:: s_n & := \frac{1}{2} (\nu - 2n - 1) \\
                  z   & := \nu \exp(-\beta x)

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
        via a more numerically stable but approximate analytic formula.

        :param x: Array of grid nodes to evaluate :math:`\mu_0` on.
        """
        # Improved version for small epsilon
        delta = ((self._nu - 1) + 1 / (12 * (self._nu - 1) - 1 / (10 * (self._nu - 1))))
        r2 = sqrt(self._beta * delta / (exp(1) * self._nu))
        r4 = ((self._nu - 1) / (2 * pi)) ** 0.25
        ex = log(delta) - 1 - log(self._nu) + (self._nu - 1) / self._nu * self._beta * x + exp(-self._beta * x)
        return r2 * r4 * exp(-self._nu / 2 * ex)


    def _grid_wrap(self, agrid):
        # TODO: Consider additional input types for "nodes":
        #       list of numpy ndarrays, list of single python scalars
        if not isinstance(agrid, AbstractGrid):
            agrid = atleast_2d(agrid)
            agrid = agrid.reshape(self._dimension, -1)
            agrid = GridWrapper(agrid)
        return agrid


    def _evaluate_basis_at_direct(self, grid):
        r"""Evaluate the eigenstates :math:`\{\mu_k\}_{k \in \mathfrak{K}}` by direct computation via the analytic formula.

        :param grid: Array of grid nodes to evaluate :math:`\mu_k` on.

        .. warning:: Evaluation becomes inaccurate of fails for large :math:`k` or
                     small :math:`\varepsilon` due to the gamma functions involved.
        """
        # The grid
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes()
        nn = grid.get_number_nodes(overall=True)

        B = zeros((self._basis_size, nn), dtype=complexfloating)

        for n in range(self._basis_size):
            B[n, :] = self._evaluate_mun(n, nodes)
        return B


    def evaluate_basis_at(self, grid):
        r"""Evaluate the eigenstates :math:`\{\mu_k\}_{k \in \mathfrak{K}}` by a recursive scheme, starting
        with an approximation to the groundstate :math:`\mu_0`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :return: An array containing the values of the :math:`\mu` at the nodes :math:`\gamma`.
        """
        # The grid
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes()
        nn = grid.get_number_nodes(overall=True)

        B = zeros((self._basis_size, nn), dtype=complexfloating)

        # Groundstate
        B[0, :] = self._evaluate_mu0(nodes)

        # Recursion
        for n in range(0, self._basis_size - 1):
            rn = sqrt((n + 1) * (self._nu - n - 1))
            ln = sqrt(n * (self._nu - n))
            sn = 1 / 2 * (self._nu - 2 * n - 1)
            pf1 =  1 / rn * sqrt((sn - 1) / sn      ) * (2 * sn)     / (2 * sn + 1)
            pf2 = ln / rn * sqrt((sn - 1) / (sn + 1)) * (2 * sn - 1) / (2 * sn + 1)
            pf1 *= ((4 * sn**2 - 1) * exp(self._beta * nodes) / self._nu - self._nu)
            # Note: Index wraps around for n = 0
            B[n + 1, :] = pf1 * B[n, :] - pf2 * B[n - 1, :]

        return B


    def evaluate_at(self, grid):
        r"""Evaluate the Morse wavepacket :math:`\mu` at the given nodes :math:`\gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :return: An array containing the values of the :math:`\mu` at the nodes :math:`\gamma`.
        """
        B = self.evaluate_basis_at(grid)
        mu = sum(self._coefficients * B, axis=0).reshape(1, -1)
        return mu


    def norm(self):
        r"""Calculate the :math:`L^2` norm :math:`\langle\mu|\mu\rangle` of the wavepacket :math:`\mu`.

        :return: The norm of :math:`\mu`.
        """
        return norm(self._coefficients)


    def energy_levels(self):
        r"""Compute the energy levels of the eigenstates :math:`\{\mu_n\}_{n \in \mathfrak{K}}` by the analytic formula:

        .. math:: E_n = -\frac{1}{2} {\left(\sqrt{2 V_0} - \varepsilon^2 |\beta| \left(n + \frac{1}{2}\right)\right)}^2
        """
        K = arange(self._basis_size)
        return -1 / 2 * (sqrt(2 * self._V0) - self._eps**2 * abs(self._beta) * (K + 1 / 2))**2
