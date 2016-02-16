"""The WaveBlocks Project

Compute the action of the gradient operator applied to
a linear combination of arbitrary wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND.Gradient import Gradient
from WaveBlocksND.LinearCombinationOfWPs import LinearCombinationOfWPs

__all__ = ["GradientLCWP"]


class GradientLCWP(Gradient):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x` applied to
    a linear combination :math:`\Upsilon` of arbitrary wavepackets.
    """

    def __init__(self):
        r"""
        """
        pass


    def apply_gradient(self, lincomb, component=None):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x`
        applied to a linear combination :math:`\Upsilon` of arbitrary wavepackets.

        :param lincomb: The lincomb :math:`\Upsilon` containing the wavepackets :math:`\Psi_j`.
        :type lincomb: A :py:class:`LinearCombinationOfWPs` instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i` of each :math:`\Psi_j`.
        :type component: Integer or ``None``.
        :return: A list of linear combinations of the gradients of all packets from the given
                 linearcombination. There are :math:`D` linearcombinations :math:`\Upsilon_d`,
                 one for each space variable component :math:`\partial_{x_d}` of the gradient.
        """
        D = lincomb.get_dimension()
        N = lincomb.get_number_components()

        # TODO: Optimizing this. For large linear combinations, computing
        #       and storing *all* gradient packets is inefficient. Maybe
        #       better go packet by packet.
        gradients = [LinearCombinationOfWPs(D, N) for d in range(D)]
        coefficients = lincomb.get_coefficients()

        for i, packet in enumerate(lincomb.get_wavepackets()):
            G = packet.get_gradient_operator()
            gradient_wps = G.apply_gradient(packet, component=component)
            for d, grad_wp in enumerate(gradient_wps):
                gradients[d].add_wavepacket(grad_wp, coefficient=coefficients[i, 0])

        return gradients
