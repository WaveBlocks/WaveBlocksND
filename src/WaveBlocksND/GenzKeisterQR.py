"""The WaveBlocks Project

This file contains the class for Genz-Keister quadrature.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
from numpy import array, indices, repeat, hstack

from QuadratureRule import QuadratureRule
from Combinatorics import partitions, lattice_points, permutations

__all__ = ["GenzKeisterQR"]


class GenzKeisterQR(QuadratureRule):
    r"""This class implements a Genz-Keister quadrature rule.
    """

    def __init__(self, dimension, level, options={}):
        r"""Initialize a new Genz-Keister quadrature rule.

        :param dimension: The dimension :math:`D` of the Genz-Keister quadrature.
        :param level: The level :math:`K` of the Genz-Keister quadrature.

        :raise: :py:class:`ValueError` if the `level` is not in the range :math:`[1, \ldots, 18]`.

        References
        ----------

        .. [1] Alan Genz and Bradley Keister: "Fully Symmetric Interpolatory Rules
               for Multiple Integrals over Infinite Regions with Gaussian Weight",
               J. Comp. Appl. Math. 71 (1996), pp. 299-309.
        """
        if not level >= 1:
            raise ValueError("Genz-Keister level has to be 1 at least.")

        if not level <= 18:
            raise ValueError("Genz-Keister quadrature rule not available for order %d" % level)

        # The space dimension of the quadrature rule.
        self._dimension = dimension

        # The level of the Genz-Keister construction.
        self._level = level

        # The order of the Genz-Keister quadrature.
        self._order = level

        # Set the options
        self._options = options

        # The number of nodes in this quadrature rule.
        self._number_nodes = None

        # The quadrature nodes \gamma.
        self._nodes = None

        # The quadrature weights \omega.
        self._weights = None


    def __str__(self):
        return "Genz-Keister quadrature rule of level %d in %d dimensions" % (self._level, self._dimension)


    def get_description(self):
        r"""Return a description of this quadrature rule object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "GenzKeisterQR"
        d["dimension"] = self._dimension
        d["order"] = self._level
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


    def _compute_nodes(self, P):
        r"""Compute fully symmetric quadrature nodes for given partition `P`.

        :param P: Partition `P`
        :return: The nodes that belong to the given partition.
        """
        D = self._dimension
        # Number of 0 entries in P
        xi = sum(P==0)
        # Compute all sign flips
        u,v = indices((D-xi, 2**(D-xi)))
        s = -2*((v>>u) & 1)+1

        nodes = []
        for P in permutations(P):
            # Generators corresponding to partition P give current node
            X = generators[P].reshape(-1,1).copy()
            # Flip signs
            I = (P != 0)
            XX = repeat(X, 2**(D-xi), 1)
            XX[I,:] *= s
            nodes.append(XX)

        return hstack(nodes)


    def _compute_weights(self, P, K):
        r"""Function to compute weights for partition :math:`P`.

        :param P: Partition :math:`P`
        :param K: Rule level :math:`K`
        :return: The weights that belong to the given partition.
        """
        D = self._dimension
        W = 0.0
        for Q in lattice_points(D, K-P.sum()):
            w = 1.0
            for i in xrange(D):
                w *= WF[P[i], Q[i]+P[i]]
            W += w

        K = sum(P>0)
        W /= 2**K
        return W


    def construct_rule(self):
        r"""Compute a Genz-Keister quadrature rule.
        """
        D = self._dimension
        K = self._level
        # Note: K is shifted by 1 wrt the original paper
        # This makes K=3 equal to the Gauss-Hermite Rule of 3 nodes.
        K -= 1

        nodes = []
        weights = []
        # Iterate over all relevant integer partitions
        for P in partitions(D, K):
            if P.sum() + Z[P].sum() <= K:
                # Compute nodes and weights for given partition
                p = self._compute_nodes(P)
                w = self._compute_weights(P, K)
                w = repeat(w, p.shape[1])
                nodes.append(p)
                weights.append(w)

        self._nodes = hstack(nodes)
        self._weights = hstack(weights)
        self._number_nodes = self._weights.size


# Magic numbers
generators = array([0.0,
                    1.2247448713915890491,
                    2.9592107790638377223,
                    0.52403354748695764515,
                    2.0232301911005156592,
                    4.4995993983103888029,
                    0.87004089535290290013,
                    3.66777421594633786,
                    1.8357079751751868738,
                    2.2665132620567880275,
                    6.3759392709822359517,
                    0.17606414208200893503,
                    5.6432578578857450628,
                    1.5794121348467670857,
                    5.0360899444730939687,
                    2.5705583765842967091,
                    4.0292201405043713648,
                    3.3491639537131949774])


WF = array([[1.7724538509055160273, -0.59081795030183867577, 0, -0.36853179227754229196, -0.36295709857235321137, 0, 0, 0, 0.022960325295871627477, 0.064774271955447564196, 0, 0, 0, 0, 0, -0.078993210637533977466, -0.21232194703927159819, -0.2456538225746588462],
            [0, 0.59081795030183867577, 0, -0.099659379601210619666, -0.15492071280527271217, 0, 0, 0, -0.021871890362310419772, -0.087151369618079604908, 0, 0, 0, 0, 0, -0.0084280964777479900247, -0.024959611548050215265, -0.033335853594384873008],
            [0, 0, 0, 0.0024661361310306921164, -0.0021319596048930450734, 0, 0, 0, 4.0299852847920501352e-05, -0.00016134458352852282233, 0, 0, 0, 0, 0, -2.2985475621124760988e-06, -1.3413229500114444752e-05, -7.0762776783582841311e-05],
            [0, 0, 0, 0.46572503574772221951, 0.49166353523456075667, 0, 0, 0, -0.054992498935111849795, -0.1639032046687536758, 0, 0, 0, 0, 0, -0.030615086786806054104, -0.083704656697661577387, -0.099275694118263682144],
            [0, 0, 0, 0, 0.028346235747958211943, 0, 0, 0, -0.0034131974823463628108, -0.047397915243750231577, 0, 0, 0, 0, 0, 0.0027516228912164496213, 0.0098895388621464664082, 0.018017220768844287104],
            [0, 0, 0, 0, 0, 0, 0, 0, 7.4926939886103445085e-08, -7.1867796345038965618e-08, 0, 0, 0, 0, 0, -2.4809313611998399183e-10, 2.6985146228426857407e-09, -3.8784905914488334667e-09],
            [0, 0, 0, 0, 0, 0, 0, 0, 0.050311651403425852408, 0.16646558769663434605, 0, 0, 0, 0, 0, 0.015121380597001994656, 0.042631804134288498327, 0.052894007283365266033],
            [0, 0, 0, 0, 0, 0, 0, 0, -2.908568677413876793e-06, 5.0691221187388293689e-06, 0, 0, 0, 0, 0, 1.7290497100670668627e-08, 2.7120025500606351173e-07, -1.5742878326131491385e-06],
            [0, 0, 0, 0, 0, 0, 0, 0, 0.0069681438693607597638, 0.057142342329530976368, 0, 0, 0, 0, 0, -0.0052916310381902192374, -0.017948738442487162366, -0.029684350935920472511],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.010226634878176754698, 0, 0, 0, 0, 0, -0.00054959323707625352453, -0.0021610363202395604912, -0.0046128892358297682251],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.4731291214089043723e-15, -7.9937156824834138633e-15, 3.5243233640534113735e-15],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.099711523786586315879, 0.26852239716792939686, 0.31153799075718904033],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5.2609392917097827721e-13, 1.4704795024993304091e-12, -9.2506568007263991711e-13],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0062420420705365663879, 0.019823703327519339518, 0.029495272419014293819],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.813505764633597492e-11, -8.6698198707035091664e-11, 7.9542508450399069479e-11],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.333031955787889199e-05, 0.00024173377847937665464, 0.00068064388805162083445],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4.7806686765607031912e-08, 1.2364713155024626226e-07],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.6925599474040720024e-06]])


Z = array([ 0,  0,  1,  0,  0,  3,  2,  1,  0,  0,  5,  4,  3,  2,  1,  0,  0,  0])
