"""The WaveBlocks Project

This file contains the abstract class for representing basis shapes.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""


class BasisShape(object):
    r"""This class defines the abstract interface to basis shapes.
    A basis shape is essentially all information and operations
    related to the set :math:`\mathcal{K}` of multi-indices :math:`k`.
    """

    def __init__(self):
        r"""
        :raise NotImplementedError: Abstract interface.
        """
        raise NotImplementedError("'BasisShape' is an abstract interface.")


    def __contains__(self, k):
        r"""
        Checks if a given multi-index :math:`k` is part of the basis set :math:`\mathcal{K}`.
        :raise NotImplementedError: Abstract interface.
        """
        raise NotImplementedError("'BasisShape' is an abstract interface.")


    def __iter__(self):
        r"""Implements iteration over the multi-indices :math:`k` of the non-extended
        basis set :math:`\mathcal{K}`.
        """
        raise NotImplementedError("'BasisShape' is an abstract interface.")


    def contains(self, k):
        r"""
        Checks if a given multi-index :math:`k` is part of the basis set :math:`\mathcal{K}`.

        :param k: The multi-index we want to test.
        :type k: tuple
        :raise NotImplementedError: Abstract interface.
        """
        return k in self


    def get_dimension(self):
        r"""
        Returns the dimension :math:`D` of the basis shape :math:`\mathcal{K}`.
        This is defined as the number of components each multi-index
        :math:`k = (k_0, \ldots, k_{D-1})` has.
        """
        return self._dimension


    def get_basissize(self, extended=False):
        r"""
        Returns the size :math:`|\mathcal{K}|` of the basis. The size is the
        number of distinct multi-indices :math:`k` that belong to the basis
        :math:`\mathcal{K}`.
        """
        if not extended:
            return self._basissize
        else:
            return self._basissize_ext


    def get_node_iterator(self, mode="lex"):
        r"""
        Returns an iterator to iterate over all basis elements :math:`k`.

        :param mode: The mode by which we iterate over the indices. Default is 'lex'
                     for lexicographical order. Supported is also 'chain', for
                     the chain-like mode, details see the manual.
        :type mode: string
        :raise NotImplementedError: Abstract interface.
        """
        raise NotImplementedError("'BasisShape' is an abstract interface.")


    def get_neighbours(self, k):
        r"""
        Returns a list of all multi-indices that are neighbours of a given
        multi-index :math:`k`. A direct neighbours is defines as
        :math:`(k_0, \ldots, k_d \pm 1, \ldots, k_{D-1}) \forall d \in [0 \ldots D-1]`.

        :param k: The multi-index of which we want to get the neighbours.
        :type k: tuple
        :return: A list containing the pairs :math:`(d, k^\prime)`.
        :raise NotImplementedError: Abstract interface.
        """
        raise NotImplementedError("'BasisShape' is an abstract interface.")
