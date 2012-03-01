r"""The WaveBlocks Project

This file contains the abstract base class for basis transformation
between the canonical basis and the basis spanned by the eigenvactors
of the potential.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import numpy

from BasisTransformation import BasisTransformation

__all__ = ["BasisTransformation"]


class BasisTransformationWF(BasisTransformation):
    """
    """

    def set_grid(self, grid):
        # Keep a reference to the grid
        self._grid = grid

        # Evaluate the eigenvectors of the potential
        self._potential.calculate_eigenvectors()
        self._ev = self._potential.evaluate_eigenvectors_at(grid)


    def transform_to_canonical(self, wavefunction):
        """Do nothing, implement a null transformation.
        """
        # No projection for potentials with a single energy level.
        # The canonical and eigenbasis are identical here.
        if self._potential.get_number_components() == 1:
            return

        # Unpack the values from the wavefunction
        values = wavefunction.get_values()

        # Reshape to match eigenvectores
        oshapes = [ value.shape for value in values ]
        values = [ value.flatten() for value in values ]

        # The basis transformation, essentially a matrix multiplication
        N = wavefunction.get_number_components()
        result = []
        for i in xrange(N):
            tmp = numpy.zeros(values[0].shape, dtype=numpy.complexfloating)
            for j in xrange(N):
                tmp += self._ev[j][i,:] * values[j]
            result.append(tmp)

        # Reshape back
        result = [ item.reshape(shape) for item, shape in zip(result, oshapes) ]

        # Pack values back to wavefunction object
        wavefunction.set_values(result)


    def transform_to_eigen(self, wavefunction):
        """Do nothing, implement a null transformation.
        """
        # No projection for potentials with a single energy level.
        # The canonical and eigenbasis are identical here.
        if self._potential.get_number_components() == 1:
            return

        # Unpack the values from the wavefunction
        values = wavefunction.get_values()

        # Reshape to match eigenvectores
        oshapes = [ value.shape for value in values ]
        values = [ value.flatten() for value in values ]

        # The basis transformation, essentially a matrix multiplication
        N = wavefunction.get_number_components()
        result = []
        for i in xrange(N):
            tmp = numpy.zeros(values[0].shape, dtype=numpy.complexfloating)
            for j in xrange(N):
                tmp += self._ev[i][j,:] * values[j]
            result.append(tmp)

        # Reshape back
        result = [ item.reshape(shape) for item, shape in zip(result, oshapes) ]

        # Pack values back to wavefunction object
        wavefunction.set_values(result)
