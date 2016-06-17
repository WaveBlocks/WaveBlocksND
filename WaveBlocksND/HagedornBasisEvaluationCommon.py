"""The WaveBlocks Project

The basic common algorithms for evaluation Hagedorn basis functions.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import atleast_2d, einsum, dot, pi
from scipy import exp
from scipy.linalg import det, inv

from WaveBlocksND.AbstractGrid import AbstractGrid
from WaveBlocksND.GridWrapper import GridWrapper

__all__ = ["HagedornBasisEvaluationCommon"]


class HagedornBasisEvaluationCommon(object):
    r"""
    """

    # We can evaluate the ground state basis function phi_0 on a set of nodes
    # the same way for homogeneous and inhomogeneous Hagedorn wavepackets.


    def _grid_wrap(self, agrid):
        # TODO: Consider additional input types for "nodes":
        #       list of numpy ndarrays, list of single python scalars
        if not isinstance(agrid, AbstractGrid):
            agrid = atleast_2d(agrid)
            agrid = agrid.reshape(self._dimension, -1)
            agrid = GridWrapper(agrid)
        return agrid


    def _evaluate_phi0(self, component, nodes, *, prefactor=False):
        r"""Evaluate the lowest order basis function :math:`\phi_0` on a
        grid :math:`\Gamma` of nodes.

        :param Pi: The parameter set :math:`\Pi`.
        :param nodes: The nodes we evaluate :math:`\phi_0` at.
        :type nodes: An ndarray of shape ``(D, |\Gamma|)``.
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :param root: The function used to compute the square root in the prefactor.
                     Defaults to the ``sqrt`` function of ``numpy`` but can be any
                     callable object and especially an instance of :py:class:`ContinuousSqrt`.
        :return: An ndarray of shape ``(|\Gamma|)``.
        """
        d = self._dimension
        eps = self._eps
        q, p, Q, P, S = self.get_parameters(component=component)

        # TODO: Use LU instead of inv(...)
        df = nodes - q
        pr1 = einsum("ik,ij,jk->k", df, dot(P, inv(Q)), df)
        pr2 = einsum("ij,ik", p, df)
        exponent = 1.0j / eps**2 * (0.5 * pr1 + pr2)

        # The problematic prefactor cancels in inner products
        if prefactor is True:
            prefactor = (pi * eps**2)**(-d * 0.25) / self._get_sqrt(component)(det(Q))
        else:
            prefactor = (pi * eps**2)**(-d * 0.25)

        return prefactor * exp(exponent)


    def evaluate_at(self, grid, *, component=None, prefactor=False):
        r"""Evaluate the Hagedorn wavepacket :math:`\Psi` at the given nodes :math:`\gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
                          (Defaults to ``None`` for evaluating all components.)
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :return: A list of arrays or a single array containing the values of the :math:`\Phi_i` at the nodes :math:`\gamma`.
        """
        Pis = self.get_parameters(component=component, aslist=True)

        if component is not None:
            phase = exp(1.0j * Pis[component][4] / self._eps**2)
            values = phase * self.slim_recursion(grid, component, prefactor=prefactor)

        else:
            values = []

            for component in range(self._number_components):
                # Note: This is very inefficient! We may evaluate the same basis functions multiple
                #       times. But as long as we don't know that the basis shapes are true subsets
                #       of the largest one, we can not evaluate just all functions in this
                #       maximal set.

                # TODO: Find more efficient way to do this

                phase = exp(1.0j * Pis[component][4] / self._eps**2)
                values.append(phase * self.slim_recursion(grid, component, prefactor=prefactor))

        return values
