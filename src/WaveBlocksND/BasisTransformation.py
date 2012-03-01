r"""The WaveBlocks Project

This file contains the abstract base class for basis transformation
between the canonical basis and the basis spanned by the eigenvactors
of the potential.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

__all__ = ["BasisTransformation"]


class BasisTransformation:
    """
    """

    def __init__(self, potential):
        """Create a new :py:class:`BasisTransformation` instance for a given potential
        matrix :math:`V(x)`.
        :param potential: The potential underlying the basis transformation.
        :type potential: A :py:class:`MatrixPotential` instance.
        """
        # Keep a reference to the potential
        self._potential = potential


    def transform_to_canonical(self, transformable):
        """Do nothing, implement a null transformation.
        """
        return transformable


    def transform_to_eigen(self, transformable):
        """Do nothing, implement a null transformation.
        """
        return transformable
