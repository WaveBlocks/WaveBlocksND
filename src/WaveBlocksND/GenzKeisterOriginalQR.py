"""The WaveBlocks Project

This file contains the class for Genz-Keister quadrature.
The quadrature is not transformed to exclude the exponential
weight factor.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
from numpy import array, indices, repeat, hstack

from QuadratureRule import QuadratureRule
from Combinatorics import partitions, lattice_points, permutations

__all__ = ["GenzKeisterOriginalQR"]


class GenzKeisterOriginalQR(QuadratureRule):
    r"""This class implements a Genz-Keister quadrature rule.
    The quadrature is not transformed to exclude the exponential
    weight factor :math:`\exp(-x^2)`.
    """

    def __init__(self, dimension, level, options={}):
        r"""Initialize a new Genz-Keister quadrature rule.

        :param dimension: The dimension :math:`D` of the Genz-Keister quadrature.
        :param level: The level :math:`K` of the Genz-Keister quadrature.

        :raise: :py:class:`ValueError` if the :math:`K` is not in the range :math:`[1, \ldots, 18]`.

        References
        ----------

        .. [1] Alan Genz: "Fully Symmetric Interpolatory Rules for Multiple Integrals",
               SIAM J. Num. Analysis. 23(6) (1986), pp. 1273-1283.

        .. [2] Alan Genz and Bradley Keister: "Fully Symmetric Interpolatory Rules
               for Multiple Integrals over Infinite Regions with Gaussian Weight",
               J. Comp. Appl. Math. 71 (1996), pp. 299-309.
        """
        # The space dimension of the quadrature rule.
        self._dimension = dimension

        # The level of the Genz-Keister construction.
        self._level = level

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
        return "Genz-Keister quadrature rule (untransformed) of level %d in %d dimensions" % (self._level, self._dimension)


    def get_description(self):
        r"""Return a description of this quadrature rule object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "GenzKeisterOriginalQR"
        d["dimension"] = self._dimension
        d["level"] = self._level
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


    def construct_rule(self, K):
        r"""Compute a Genz-Keister quadrature rule.

        :param K: The level :math:`K` of the Genz-Keister construction.

        .. note:: This is an internal method and there should be no reason
                  to explicitely call it manually.

        .. warning:: This method can be expensive and may take some time
                     to finish. Also, the quadrature nodes may use large amounts
                     of memory depending on the dimension and level parameters.
        """
        # Check for valid level parameter
        if not K >= 1:
            raise ValueError("Genz-Keister level has to be 1 at least.")
        if not K <= 30:
            raise ValueError("Genz-Keister quadrature rule not available for level %d" % K)

        self._level = K
        D = self._dimension

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

        self._nodes = hstack(nodes).reshape(D,-1)
        self._weights = hstack(weights).reshape(1,-1)
        self._number_nodes = self._weights.size


# Magic numbers
# For derivation and computation see our technical report.
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
                    3.3491639537131949774,
                    12.371183263294440156,
                    0.36668252574926773363,
                    11.773315693849850411,
                    0.66761453794663251987,
                    11.279571841264790728,
                    1.0853772883690724485,
                    10.839884501585234819,
                    1.3554874833640409297,
                    10.435144794449726187,
                    1.8804002593778771426,
                    10.055514590896118546,
                    2.4894835291142853745])

WF = array([[1.7724538509055160273, -0.59081795030183867577, 0, -0.36853179227754229196, -0.36295709857235321137, 0, 0, 0, 0.022960325295871627477, 0.064774271955447564196, 0, 0, 0, 0, 0, -0.078993210637533977466, -0.21232194703927159819, -0.2456538225746588462, 0, 0, 0, 0, 0, 0, 0, 0, 0.0069854780851807704795, 0.93317381243987030624, 1.4408446595539685687, 32.462751316915930367],
            [0, 0.59081795030183867577, 0, -0.099659379601210619666, -0.15492071280527271217, 0, 0, 0, -0.021871890362310419772, -0.087151369618079604908, 0, 0, 0, 0, 0, -0.0084280964777479900247, -0.024959611548050215265, -0.033335853594384873008, 0, 0, 0, 0, 0, 0, 0, 0, -0.00083442479578140975084, -0.19359613215029009931, -0.30341863259901277748, -9.0190236387036132044],
            [0, 0, 0, 0.0024661361310306921164, -0.0021319596048930450734, 0, 0, 0, 4.0299852847920501352e-05, -0.00016134458352852282233, 0, 0, 0, 0, 0, -2.2985475621124760988e-06, -1.3413229500114444752e-05, -7.0762776783582841311e-05, 0, 0, 0, 0, 0, 0, 0, 0, 9.9226143816067400596e-11, -8.9771296128557796343e-09, -1.5175167948853198633e-08, 8.2790782248066394335e-07],
            [0, 0, 0, 0.46572503574772221951, 0.49166353523456075667, 0, 0, 0, -0.054992498935111849795, -0.1639032046687536758, 0, 0, 0, 0, 0, -0.030615086786806054104, -0.083704656697661577387, -0.099275694118263682144, 0, 0, 0, 0, 0, 0, 0, 0, -0.010933454407877926348, -1.583559780719271216, -2.4517160773325681828, -57.79911654719117283],
            [0, 0, 0, 0, 0.028346235747958211943, 0, 0, 0, -0.0034131974823463628108, -0.047397915243750231577, 0, 0, 0, 0, 0, 0.0027516228912164496213, 0.0098895388621464664082, 0.018017220768844287104, 0, 0, 0, 0, 0, 0, 0, 0, -8.2396289613717090217e-07, 0.00069805033976115270287, 0.0011232827256560501634, 0.074544571196262348206],
            [0, 0, 0, 0, 0, 0, 0, 0, 7.4926939886103445085e-08, -7.1867796345038965618e-08, 0, 0, 0, 0, 0, -2.4809313611998399183e-10, 2.6985146228426857407e-09, -3.8784905914488334667e-09, 0, 0, 0, 0, 0, 0, 0, 0, 2.4377187921244046447e-16, -6.8906698719329670782e-15, -1.3303118025999139614e-14, 1.3222051263169226208e-13],
            [0, 0, 0, 0, 0, 0, 0, 0, 0.050311651403425852408, 0.16646558769663434605, 0, 0, 0, 0, 0, 0.015121380597001994656, 0.042631804134288498327, 0.052894007283365266033, 0, 0, 0, 0, 0, 0, 0, 0, -0.00228042374526934975, -0.38761835932549071167, -0.60300716658807476412, -15.476251566409639582],
            [0, 0, 0, 0, 0, 0, 0, 0, -2.908568677413876793e-06, 5.0691221187388293689e-06, 0, 0, 0, 0, 0, 1.7290497100670668627e-08, 2.7120025500606351173e-07, -1.5742878326131491385e-06, 0, 0, 0, 0, 0, 0, 0, 0, 4.1133482077574836997e-13, -1.9592820190556365775e-11, -3.4894317864083047118e-11, 6.7158596618593598195e-10],
            [0, 0, 0, 0, 0, 0, 0, 0, 0.0069681438693607597638, 0.057142342329530976368, 0, 0, 0, 0, 0, -0.0052916310381902192374, -0.017948738442487162366, -0.029684350935920472511, 0, 0, 0, 0, 0, 0, 0, 0, 3.9407676130575873034e-06, 0.01120796927808146534, 0.017902020219827321149, 0.88400496224987437599],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.010226634878176754698, 0, 0, 0, 0, 0, -0.00054959323707625352453, -0.0021610363202395604912, -0.0046128892358297682251, 0, 0, 0, 0, 0, 0, 0, 0, 6.8201331441269281092e-08, -2.0119651796767681865e-05, -3.272801877689525844e-05, -0.0043094107121951234471],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.4731291214089043723e-15, -7.9937156824834138633e-15, 3.5243233640534113735e-15, 0, 0, 0, 0, 0, 0, 0, 0, -3.5979073644975835664e-23, 4.5787578394467178907e-22, 1.1823259066197271028e-21, -4.7914913139336527046e-21],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.099711523786586315879, 0.26852239716792939686, 0.31153799075718904033, 0, 0, 0, 0, 0, 0, 0, 0, -0.012942205205303980778, -1.7442103174041754507, -2.6939320223251190934, -61.000372377628353199],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5.2609392917097827721e-13, 1.4704795024993304091e-12, -9.2506568007263991711e-13, 0, 0, 0, 0, 0, 0, 0, 0, 1.5785932351677890194e-20, -2.6338462204337911645e-19, -5.9364547960338525453e-19, 3.2318133244861044316e-18],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0062420420705365663879, 0.019823703327519339518, 0.029495272419014293819, 0, 0, 0, 0, 0, 0, 0, 0, -2.8712227378685372603e-05, -0.013023628116982740502, -0.020617469269662274043, -0.77744561213595182479],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.813505764633597492e-11, -8.6698198707035091664e-11, 7.9542508450399069479e-11, 0, 0, 0, 0, 0, 0, 0, 0, -2.4991256866116674136e-18, 5.4084778438692794698e-17, 1.1146766899167308562e-16, -8.1214610992049753868e-16],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.333031955787889199e-05, 0.00024173377847937665464, 0.00068064388805162083445, 0, 0, 0, 0, 0, 0, 0, 0, -3.1728953354017843925e-09, 4.8788852162859657901e-07, 8.0598370444025131215e-07, -0.00027432960967869342032],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4.7806686765607031912e-08, 1.2364713155024626226e-07, 0, 0, 0, 0, 0, 0, 0, 0, -1.6391871509412773069e-14, 6.0972774207556211602e-13, 1.1215022650753154351e-12, -1.5601959912306827249e-11],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.6925599474040720024e-06, 0, 0, 0, 0, 0, 0, 0, 0, -5.0478028338888555907e-12, 3.1042162446656135191e-10, 5.3910396570270701411e-10, -1.4997190640433153694e-08],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.966610468031952332e-35, -9.3725229288180246934e-35, 2.8175911769125794736e-34, -2.6791356729768828299e-34],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.012234307825621222445, 1.6989572266556452903, 2.6267270916601100203, 60.493525291298403679],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.7462006007817010057e-33, 6.1064103564975326083e-33, -2.5424116143210707669e-32, 2.6810227097986687364e-32],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0062338782069609262528, 0.95288239924914269787, 1.4777893468529787418, 35.875172109135358315],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.6404631535319264035e-32, -1.0083297438305276761e-31, 6.0279502653589950369e-31, -6.9543885865056008569e-31],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0012910899749692758878, 0.25864570545974318927, 0.40406333226476557602, 11.240276946785162355],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.4853386161233029025e-31, 6.1561988785958881913e-31, -5.8641391245757267453e-30, 7.35654767983557935e-30],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00027128436115260543027, 0.075441650507782930325, 0.11863969202649792811, 3.7993742821956509655],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.8439821915492598744e-31, -1.2750675795119443774e-30, 2.5590556531634099977e-29, -3.4795116876601401699e-29],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0089789557648437572205, -0.014366120484444668385, -0.75366721653557697842],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4.2040208044811380297e-29, 6.18461706215374787e-29],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00081040558279177164563]])

Z =  array([ 0, 0, 1, 0, 0, 3, 2, 1, 0, 0, 5, 4, 3, 2, 1, 0, 0, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0])
