"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import add

from WaveBlocksND.WaveFunction import WaveFunction
from WaveBlocksND.BlockFactory import BlockFactory
from functools import reduce

__all__ = ["Initializer"]


class Initializer(object):

    def __init__(self, parameters):
        self._parameters = parameters


    def initialize_for_fourier(self, grid):
        # NOTE: This code is only temporary until we get the
        #       more versatile initial value specifications.

        # TODO: Make a specification of the IV setup in configurations
        Psi = []

        for packet_descr in self._parameters["initvals"]:
            packet = BlockFactory().create_wavepacket(packet_descr)

            # Evaluate the
            X = grid.get_nodes(flat=True)
            values = packet.evaluate_at(X, prefactor=True)

            # Reshape values into hypercubic shape
            values = [ val.reshape(grid.get_number_nodes()) for val in values ]
            Psi.append(values)

        # TODO: Maybe sum up immediately instead of at the end to reduce memory usage
        Psi = reduce(add, Psi)

        # Pack the values in a WaveFunction instance
        WF = WaveFunction(self._parameters)
        WF.set_grid(grid)
        WF.set_values(Psi)

        return WF
