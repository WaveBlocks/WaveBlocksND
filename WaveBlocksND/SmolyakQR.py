"""The WaveBlocks Project

This file contains the class for constructing high dimensional
quadrature rules from one-dimensional ones by using the
Smolyak construction.

@author: R. Bourquin
@copyright: Copyright (C) 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
import operator as op
from functools import reduce
from numpy import hstack, vsplit, squeeze, lexsort, where, vstack, multiply
from numpy.linalg import norm
from scipy.special import binom

from WaveBlocksND.QuadratureRule import QuadratureRule
from WaveBlocksND.Combinatorics import lattice_points_norm
from WaveBlocksND.Utils import meshgrid_nd


__all__ = ["SmolyakQR"]


class SmolyakQR(QuadratureRule):
    r"""This class implements the construction of high dimensional
    quadrature rules from one-dimensional ones by applying the
    Smolyak construction.
    """

    def __init__(self, dimension, level, rules, options={}):
        r"""
        Initialize a :py:class:`SmolyakQR` instance.

        :param dimension: The dimension :math:`D` of the final
                          quadrature rule.
        :param level: The level :math:`k` of the Smolyak construction.
                      From theory we know that a Smolyak rule of order
                      :math:`k` is exact up to :math:`2 k - 1` if the
                      individual rules :math:`Q_i` are exact up to
                      :math:`2 i - 1`. The level has to be larger than
                      or equal 1.
        :param rules: A collection of :py:class:`QuadratureRule` subclass
                      instances. Their nodes and weights will be used
                      in the Smolyak construction.

        .. note:: The ``rules`` object must implement simple indexing by
                  non-negative numbers and for a given index :math:`j` return
                  a univariate quadrature rule of order :math:`j`. Special attention
                  must be payed in case this object is mutable. We do not copy it.

        .. warning:: This implementation uses a special optimization that
                     speeds up the computation of all quadrature nodes, especially
                     in high dimension, but is only valid if all one dimensional
                     rules have symmetric nodes and weights. That is, for every node
                     weight pair :math:`(\gamma, \omega)` that is part of the rule,
                     the pair :math:`(-\gamma, \omega)` is also contained
                     in the quadrature rule.
        """
        # The space dimension of the quadrature rule.
        self._dimension = dimension

        # The level of the Smolyak construction.
        self._level = level

        # The individual quadrature rules.
        self._rules = rules

        # Set the options
        self._options = options

        # The number of nodes in this quadrature rule.
        self._number_nodes = None

        # The quadrature nodes \gamma.
        self._nodes = None

        # The quadrature weights \omega.
        self._weights = None

        # Actually compute the nodes and weights.
        self.construct_rule(level)


    def __str__(self):
        s = "Sparse grid (Smolyak) quadrature rule consisting of:\n"
        l = ["  " + str(rule) + "\n" for k,rule in self._rules.items() if k <= self._level]
        s += reduce(op.add, l)
        return s


    def get_description(self):
        r"""Return a description of this quadrature rule object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "SmolyakQR"
        d["dimension"] = self._dimension
        d["qr_rules"] = [ qr.get_description() for qr in self._rules ]
        d["options"] = deepcopy(self._options)
        return d


    def get_nodes(self, flat=True, split=False):
        r"""Returns the quadrature nodes :math:`\{\gamma_i\}_i`.

        :param flat: Dummy parameter for API compatibility with Grids.
        :type flat: Boolean, default is ``True``.
        :param split: Dummy parameter for API compatibility with Grids.
        :type split: Boolean, default is ``False``.
        :return: An ndarray containing the quadrature nodes :math:`\{\gamma_i\}_i`.
        """
        if self._nodes is None:
            self.construct_rule()
        return self._nodes.copy()


    def get_weights(self):
        r"""Returns the quadrature weights :math:`\{\omega_i\}_i`.

        :return: An ndarray containing the quadrature weights :math:`\{\omega_i\}_i`.
        """
        if self._weights is None:
            self.construct_rule()
        return self._weights.copy()


    def construct_rule(self, K, tolerance=1e-15):
        r"""Compute the quadrature nodes :math:`\{\gamma_i\}_i` and quadrature
        weights :math:`\{\omega_i\}_i`.

        :param K: The level :math:`K` of the Smolyak construction.
        :param tolerance: Tolerance for dropping identical quadrature nodes.

        .. note:: This is an internal method and there should be no reason
                  to explicitely call it manually.

        .. warning:: This method is very expensive and may take a long time
                     to finish. Also, the quadrature nodes may use large amounts
                     of memory depending on the dimension and level parameters.
        """
        # Check for valid level parameter
        if not K >= 1:
            raise ValueError("Smolyak level has to be 1 at least.")
        if K > max(self._rules.keys()):
            raise ValueError("Not enough quadrature rules to build Smolyak grid of level %d" % K)

        self._level = K
        D = self._dimension

        allnodes = []
        allweights = []
        factors = []

        # Index Set
        for q in range(max(0, K-D), K):
            S = lattice_points_norm(D, q)
            for j, s in enumerate(S):
                # Only use non-negative nodes for the construction.
                # The quadrature nodes \gamma.
                rules = [ self._rules[si+1] for si in s ]
                nodes = [ rule.get_nodes() for rule in rules ]
                indices = [ where(n >= 0) for n in nodes ]
                nodes = meshgrid_nd([ n[i] for n, i in zip(nodes, indices) ])
                nodes = vstack([ node.reshape(-1) for node in nodes ])
                # The quadrature weights \omega.
                weights = [ rule.get_weights() for rule in rules ]
                weights = meshgrid_nd([ w[i] for w, i in zip(weights, indices) ])
                weights = reduce(multiply, weights)
                allnodes.append(nodes)
                allweights.append(weights.reshape(-1))
                factors.append((-1)**(K-1-q) * binom(D-1, K-1-q))

        # Sort
        allnodes = hstack(allnodes).reshape(D,-1)
        allweights = hstack([f*w for f, w in zip(factors,allweights)]).reshape(1,-1)

        I = lexsort(map(squeeze, vsplit(allnodes, D))[::-1])

        allnodes = allnodes[:,I].reshape(D,-1)
        allweights = allweights[:,I].reshape(1,-1)

        # Remove duplicates
        last = 0
        I = [last]

        no = norm(allnodes[:,:-1] - allnodes[:,1:], axis=0)

        for col in range(1, allnodes.shape[1]):
            if no[col-1] <= tolerance:
                allweights[0,last] += allweights[0,col]
                allweights[0,col] = 0
            else:
                last = col
                I.append(last)

        allnodes = allnodes[:,I]
        allweights = allweights[:,I]

        # Mirror points to all other hyperoctants
        for d in range(D):
            indices = abs(allnodes[d,:]) >= tolerance
            mirrorn = allnodes[:,indices]
            mirrorn[d,:] *= -1.0
            mirrorw = allweights[:,indices]
            allnodes = hstack([allnodes, mirrorn]).reshape(D,-1)
            allweights = hstack([allweights, mirrorw]).reshape(1,-1)

        self._nodes = allnodes
        self._weights = allweights
        self._number_nodes = allnodes.shape[1]
