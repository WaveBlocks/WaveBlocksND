"""The WaveBlocks Project

This file contains code for the representation of the kinetic
operator :math:`-\frac{1}{2} \Delta` in the Fourier space.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import roll, exp
from numpy.fft import fftfreq

__all__ = ["KineticOperator"]


class KineticOperator(object):
    """
    """

    def __init__(self, grid):
        # TODO: Consider class "FourierSpaceGrid"
        D = grid.get_dimension()
        N = grid.get_number_nodes()

        shape = [-1] + (D-1)*[1]

        omega = [ fftfreq(n, d=1.0/n) for n in N ]
        omega = [ o.reshape(roll(shape, i)) for i, o in enumerate(omega) ]
        omega_sqr = map(lambda x: x*x, omega)
        omega_sqr = reduce(lambda x,y: x+y, omega_sqr)

        # Fourier space grid axes
        self._omega = omega
        # omega dot omega
        self._omega_sqr = omega_sqr


    def get_fourier_grid_axes(self):
        return tuple([ o.copy() for o in self._omega ])


    # TODO: Reconsider external API: which version of the
    #       calculate function do we want to use?
    #       Consistency with the potential?

    def calculate_operator(self, dt, eps):
        self._dt = dt
        self._eps = eps
        self._OpT = 0.5 * self._eps**4 * self._omega_sqr
        self._expOpT = exp(-0.5j * self._dt * self._eps**2 * self._omega_sqr)


    def calculate_exponential(self, factor=1):
        #self._expOpT = numpy.exp(-0.5j * PA["dt"] * PA["eps"]**2 * O)
        # TODO: Should we include the -0.5j ?
        pass


    def evaluate_at(self, grid=None):
        return self._OpT.copy()


    def evaluate_exponential_at(self, grid=None):
        return self._expOpT.copy()
