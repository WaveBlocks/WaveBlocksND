"""The WaveBlocks Project

This file contains the class for Genz-Keister quadrature.

@author: R. Bourquin
@copyright: Copyright (C) 2014, 2015 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy
from numpy import exp
from numpy.linalg import norm

from WaveBlocksND.GenzKeisterOriginalQR import GenzKeisterOriginalQR

__all__ = ["GenzKeisterQR"]


class GenzKeisterQR(GenzKeisterOriginalQR):
    r"""This class implements a Genz-Keister quadrature rule.
    """

    def __init__(self, dimension, level, options={}):
        r"""Initialize a new Genz-Keister quadrature rule.

        :param dimension: The dimension :math:`D` of the Genz-Keister quadrature.
        :param level: The level :math:`K` of the Genz-Keister quadrature.

        :raise: :py:class:`ValueError` if the :math:`K` is not in the range :math:`[1, \ldots, 18]`.

        .. [1] Alan Genz: "Fully Symmetric Interpolatory Rules for Multiple Integrals",
               SIAM J. Num. Analysis. 23(6) (1986), pp. 1273-1283.

        .. [2] Alan Genz and Bradley Keister: "Fully Symmetric Interpolatory Rules
               for Multiple Integrals over Infinite Regions with Gaussian Weight",
               J. Comp. Appl. Math. 71 (1996), pp. 299-309.
        """
        GenzKeisterOriginalQR.__init__(self, dimension, level, options=options)
        # Transform weights
        # TODO: This is the best transform we can do right now
        self._weights /= exp(-norm(self._nodes, axis=0)**2)


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
        d["level"] = self._level
        d["options"] = deepcopy(self._options)
        return d
