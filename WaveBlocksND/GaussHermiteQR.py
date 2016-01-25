"""The WaveBlocks Project

This file contains the class for Gauss-Hermite quadrature.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
from orthogonal import psi_roots

from QuadratureRule import QuadratureRule

__all__ = ["GaussHermiteQR"]


class GaussHermiteQR(QuadratureRule):
    r"""This class implements a Gauss-Hermite quadrature rule.
    """

    def __init__(self, order, options={}):
        r"""Initialize a new quadrature rule.

        :param order: The order :math:`k` of the Gauss-Hermite quadrature.
                      From theory we know that a Gauss quadrature rule
                      of order :math:`k` is exact for polynomials up to
                      degree :math:`2 k - 1`.

        :raise: :py:class:`ValueError` if the ``order`` is not 1 or above.
        """
        # Quadrature has to have at least a single (node,weight) pair.
        if not order > 0:
            raise ValueError("Quadrature rule has to be of order 1 at least.")

        # The space dimension of the quadrature rule.
        self._dimension = 1

        # The order of the Gauss-Hermite quadrature.
        self._order = order

        # Set the options
        self._options = options

        nodes, weights = psi_roots(self._order)

        # The number of nodes in this quadrature rule
        self._number_nodes = nodes.size

        # The quadrature nodes \gamma.
        self._nodes = nodes.reshape((1,self._number_nodes))

        # The quadrature weights \omega.
        self._weights = weights.reshape((1,self._number_nodes))


    def __str__(self):
        return "Gauss-Hermite quadrature rule of order %d" % self._order


    def get_description(self):
        r"""Return a description of this quadrature rule object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "GaussHermiteQR"
        d["dimension"] = self._dimension
        d["order"] = self._order
        d["options"] = deepcopy(self._options)
        return d


    def get_nodes(self):
        r"""Returns the quadrature nodes :math:`\{\gamma_i\}_i`.

        :return: An array containing the quadrature nodes :math:`\{\gamma_i\}_i`.
        """
        return self._nodes.copy()


    def get_weights(self):
        r"""Returns the quadrature weights :math:`\{\omega_i\}_i`.

        :return: An array containing the quadrature weights :math:`\{\omega_i\}_i`.
        """
        return self._weights.copy()
