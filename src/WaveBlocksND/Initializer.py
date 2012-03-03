"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import array, zeros

from WaveFunction import WaveFunction


class Initializer(object):

    def __init__(self, parameters):
        self._parameters = parameters


    def initialize_for_fourier(self, grid):
        # NOTE: This code is only temporary until we get the
        #       more versatile initial value specifications.

        # TODO: Make a specification of the IV setup in configurations

        # Import a phi_{0,0,...} wavepacket
        from phi0 import hawp

        # Set up values for a one-component wavefunction
        dim = self._parameters["dimension"]

        Pi = [ self._parameters["q"],
               self._parameters["p"],
               self._parameters["Q"],
               self._parameters["P"] ]
        Pi = map(array, Pi)

        eps = self._parameters["eps"]

        X = grid.get_nodes()

        # Evaluate the phi_0
        values = hawp(dim, Pi, eps, X)

        # Reshape values into hypercubic shape
        values = values.reshape(grid.get_number_nodes())

        psi = [ values, zeros(values.shape) ]

        # Pack the values in a WaveFunction instance
        WF = WaveFunction(self._parameters)
        WF.set_grid(grid)
        WF.set_values(psi)

        return WF
