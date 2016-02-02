"""The WaveBlocks Project

This file contains the class for representing the simplex
basis shape which is a special type of sparse basis set.

@author: R. Bourquin
@copyright: Copyright (C) 2015 R. Bourquin
@license: Modified BSD License
"""

from numpy import eye, vstack, integer

from .BasisShape import BasisShape

__all__ = ["SimplexShape"]


class SimplexShape(BasisShape):
    r"""This class implements the simplex basis shape which is a special type
    of sparse basis set. A basis shape is essentially all information and operations
    related to the set :math:`\mathfrak{K}` of multi-indices :math:`k`. The simplex
    shape in :math:`D` dimensions and with maximal 1-norm :math:`K` is defined as the set

    .. math::
        \mathfrak{K}(D, K) := \{ (k_0, \ldots, k_{D-1}) |
                                 \sum_{d=0}^{D-1} k_d = \| k \|_1 \leq K \}
    """

    def __init__(self, D, K):
        r"""
        :param D: The dimension :math:`D`
        :param K: The maximal 1-norm :math:`K`
        """
        # The dimension of K
        self._dimension = D

        # The maximal 1-norm parameter
        self._maxnorm = K

        # The linear mapping k -> index for the basis
        iil = self._get_index_iterator_lex()
        self._lima = {k:index for index, k in enumerate(iil)}
        # And the inverse mapping
        self._lima_inv = {v:k for k, v in self._lima.items()}

        # The basis size
        self._basissize = len(self._lima)


    def __str__(self):
        r""":return: A string describing the basis shape :math:`\mathfrak{K}`.
        """
        s = ("Simplex basis shape of dimension "+str(self._dimension)+
             " and maximal 1-norm "+str(self._maxnorm)+".")
        return s


    def __hash__(self):
        r"""Compute a unique hash for the basis shape. In the case of simplex
        basis shapes :math:`\mathfrak{K}` the basis is fully specified by its
        dimension :math:`D` and the maximal 1-norm parameter :math:`K`.
        """
        return hash(("SimplexShape", self._dimension, self._maxnorm))


    def __getitem__(self, k):
        r"""Make map look ups.
        """
        if type(k) is tuple or type(k) is list:
            k = tuple(k)
            assert len(k) == self._dimension
            if k in self._lima:
                return self._lima[k]
        elif type(k) is int:
            if k in self._lima_inv:
                return self._lima_inv[k]
        else:
            raise IndexError("Wrong index type")


    def __contains__(self, k):
        r"""
        Checks if a given multi-index :math:`k` is part of the basis set :math:`\mathfrak{K}`.

        :param k: The multi-index :math:`k` we want to test.
        :type k: tuple
        """
        assert len(tuple(k)) == self._dimension
        return tuple(k) in self._lima


    def __iter__(self):
        r"""Implements iteration over the multi-indices :math:`k`
        of the basis set :math:`\mathfrak{K}`.

        Note: The order of iteration is NOT fixed. If you need a special
        iteration scheme, use :py:meth:`get_node_iterator`.
        """
        # TODO: Better remove this as it may cause unexpected behaviour?
        return iter(self._lima)


    def contains(self, k):
        r"""
        Checks if a given multi-index :math:`k` is part of the basis set :math:`\mathfrak{K}`.

        :param k: The multi-index :math:`k` we want to test.
        :type k: tuple
        """
        return tuple(k) in self._lima


    def get_description(self):
        r"""Return a description of this basis shape object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current basis shape. A description
        never contains any data.
        """
        d = {}
        d["type"] = "SimplexShape"
        d["dimension"] = self._dimension
        d["K"] = self._maxnorm
        return d


    def extend(self):
        r"""Extend the basis shape such that (at least) all neighbours of all
        boundary nodes are included in the extended basis shape.
        """
        D = self._dimension
        K = self._maxnorm
        return SimplexShape(D, K+1)


    def _get_index_iterator_lex(self):
        r"""
        """
        # The maximal 1-norm
        Kmax = self._maxnorm

        def index_iterator_lex(Kmax):
            # Initialize a counter
            z = [0 for i in range(self._dimension + 1)]

            while z[self._dimension] == 0:
                # Yield the current index vector
                yield tuple(reversed(z[:-1]))

                # Increment fastest varying bit
                z[0] += 1

                # Reset overflows
                for d in range(self._dimension):
                    K = sum(z[:-1])
                    if K > Kmax:
                        z[d] = 0
                        z[d+1] += 1

        return index_iterator_lex(Kmax)


    def _get_index_iterator_chain(self, direction=0):
        r"""
        """
        def index_iterator_chain(Kmax, d):
            D = self._dimension
            # The counter
            z = [ 0 for i in range(D + 1) ]

            # Iterate over all valid stencil points
            while z[D] == 0:
                yield tuple(reversed(z[:-1]))

                # Increase index in the dimension we build the chain
                z[D-d-1] += 1

                # Check if we are done with the current base point
                # If yes, move base point and start a new chain
                # Reset overflows
                for i in range(D-d-1, D):
                    K = sum(z[(D-d-1):-1])
                    if K > Kmax:
                        z[i] = 0
                        z[i+1] += 1

        return index_iterator_chain(self._maxnorm, direction)


    def _get_index_iterator_mag(self):
        r"""
        """
        # Nodes sorted by l_1 magnitude
        nodes = sorted(self._lima.keys(), key=sum)

        def index_iterator_mag(nodes):
            for node in nodes:
                yield node

        return index_iterator_mag(nodes)


    def get_node_iterator(self, mode="lex", direction=None):
        r"""
        Returns an iterator to iterate over all basis elements :math:`k \in \mathfrak{K}`.

        :param mode: The mode by which we iterate over the indices. Default is ``lex``
                     for lexicographical order. Supported is also ``chain``, for
                     the chain-like mode, details see the manual.
        :type mode: string
        :param direction: If iterating in `chainmode` this specifies the direction
                          the chains go.
        :type direction: integer.
        """
        if mode == "lex":
            return self._get_index_iterator_lex()
        elif mode == "chain":
            if direction < self._dimension:
                return self._get_index_iterator_chain(direction=direction)
            else:
                raise ValueError("Can not build iterator for this direction.")
        elif mode == "mag":
            return self._get_index_iterator_mag()
        # TODO: Consider boundary node only iterator
        else:
            raise ValueError("Unknown iterator mode: "+str(mode)+".")


    def get_limits(self):
        r"""Returns the upper limit :math:`K` which is the same for all directions :math:`d`.

        :return: A tuple of the maximum of the multi-index in each direction.
        """
        return tuple(self._dimension * [self._maxnorm])


    def get_neighbours(self, k, selection=None, direction=None):
        r"""
        Returns a list of all multi-indices that are neighbours of a given
        multi-index :math:`k`. A direct neighbour is defined as
        :math:`(k_0, \ldots, k_d \pm 1, \ldots, k_{D-1}) \forall d \in [0 \ldots D-1]`.

        :param k: The multi-index of which we want to get the neighbours.
        :type k: tuple
        :param selection:
        :type selection: string with fixed values ``forward``, ``backward`` or ``all``.
                         The values ``all`` is equivalent to the value ``None`` (default).
        :param direction: The direction :math:`0 \leq d < D` in which we want to find
                          the neighbours :math:`k \pm e_d`.
        :type direction: int
        :return: A list containing the pairs :math:`(d, k^\prime)`.
        """
        assert len(tuple(k)) == self._dimension

        # First build a list of potential neighbours
        I = eye(self._dimension, dtype=integer)
        ki = vstack(k)

        # Forward and backward direct neighbours
        nbfw = ki + I
        nbbw = ki - I

        # Keep only the valid ones
        nbh = []

        if direction is not None:
            directions = [ direction ]
        else:
            directions = range(self._dimension)

        for d in directions:
            nfw = tuple(nbfw[:,d])
            nbw = tuple(nbbw[:,d])

            # TODO: Try to simplify these nested if blocks
            if selection in ("backward", "all", None):
                if not k[d] == 0:
                    nbh.append((d, nbw))

            if selection in ("forward", "all", None):
                if not sum(k) == self._maxnorm:
                    nbh.append((d, nfw))

        return nbh
