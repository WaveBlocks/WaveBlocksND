"""The WaveBlocks Project

This file contains the class for constructing high dimensional
quadrature rules from one-dimensional ones by taking tensor products.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
import operator as op
from functools import reduce
from numpy import vstack, multiply

from WaveBlocksND.QuadratureRule import QuadratureRule
from WaveBlocksND.Utils import meshgrid_nd

__all__ = ["TensorProductQR"]


class TensorProductQR(QuadratureRule):
    r"""This class implements the construction of high dimensional
    quadrature rules from one-dimensional ones by taking tensor products.
    """

    def __init__(self, rules, options={}):
        r"""
        Initialize a :py:class:`TensorProductQR` instance.

        :param rules: A list of :py:class:`QuadratureRule` subclass instances. Their
                      nodes and weights will be used to compute the tensor product.
        """
        if len(rules) < 1:
            raise ValueError("Can not build tensor product rule with no rules given.")

        # The dimension of the quadrature rule.
        self._dimension = len(rules)

        # The individual quadrature rules.
        self._rules = tuple(rules)

        # Set the options
        self._options = options

        # The number of nodes in this quadrature rule.
        self._number_nodes = None

        # The quadrature nodes \gamma.
        self._nodes = None

        # The quadrature weights \omega.
        self._weights = None

        # Actually compute the nodes and weights.
        self.construct_rule()


    def construct_rule(self):
        r"""Compute the tensor product of the given quadrature rules.

        .. note:: This is an internal method and there should be no reason
                  to explicitely call it manually.

        :return: The nodes :math:`\{\gamma_i\}_i` and weights :math:`\{\omega_i\}_i`
                 of the tensor product quadrature rule. The array of all
                 nodes has a shape of :math:`(D, |\Gamma|)` and the
                 array of weights is of shape :math:`(|\Gamma|)`.
        """
        self._number_nodes = reduce(multiply, [ rule.get_number_nodes() for rule in self._rules ])
        # The quadrature nodes \gamma.
        nodes = meshgrid_nd([ rule.get_nodes() for rule in self._rules ])
        self._nodes = vstack([ node.reshape(1,-1) for node in nodes ])
        # The quadrature weights \omega.
        weights = meshgrid_nd([ rule.get_weights() for rule in self._rules ])
        weights = reduce(multiply, weights)
        self._weights = weights.reshape(1,-1)


    def __str__(self):
        s = "Tensor product quadrature rule consisting of:\n"
        l = ["  " + str(rule) + "\n" for rule in self._rules]
        s += reduce(op.add, l)
        return s


    def get_description(self):
        r"""Return a description of this quadrature rule object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "TensorProductQR"
        d["dimension"] = self._dimension
        d["qr_rules"] = [ qr.get_description() for qr in self._rules ]
        d["options"] = deepcopy(self._options)
        return d


    def get_nodes(self, flat=True, split=False):
        r"""Return the quadrature nodes :math:`\{\gamma_i\}_i`.

        :param flat: Dummy parameter for API compatibility with Grids.
        :type flat: Boolean, default is ``True``.
        :param split: Dummy parameter for API compatibility with Grids.
        :type split: Boolean, default is ``False``.
        :return: An :py:class:`ndarray` containing the quadrature nodes :math:`\{\gamma_i\}_i`.
        """
        return self._nodes.copy()


    def get_weights(self):
        r"""Returns the quadrature weights :math:`\{\omega_i\}_i`.

        :return: An :py:class:`ndarray` containing the quadrature weights :math:`\{\omega_i\}_i`.
        """
        return self._weights.copy()
