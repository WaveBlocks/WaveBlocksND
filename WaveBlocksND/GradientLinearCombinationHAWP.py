"""The WaveBlocks Project

Compute the action of the gradient operator applied to a
linear combination of Hagedorn wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import squeeze

from .Gradient import Gradient
from .GradientHAWP import GradientHAWP
from .LinearCombinationOfHAWPs import LinearCombinationOfHAWPs

__all__ = ["GradientLinearCombinationHAWP"]


class GradientLinearCombinationHAWP(Gradient):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x` applied to
    a linear combination :math:`\Upsilon` of Hagedorn wavepackets :math:`\Psi`.
    """

    def __init__(self):
        r"""
        """
        pass


    # TODO: Find a more efficient way to compute gradients

    def apply_gradient(self, lincomb, component=None):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x`
        on the linear combination :math:`\Upsilon` of Hagedorn wavepackets :math:`\Psi`.

        :param lincomb: The linear combination :math:`\Upsilon`.
        :type lincomb: A :py:class:`LinearCombinationOfHAWPs` instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer or ``None``.
        :return: One linear combination :math:`\Upsilon_d` containing the gradients
                 for the component :math:`\partial_{x_d}` for each space dimension
                 component :math:`d = 1, \ldots, D`.
        """
        D = lincomb.get_dimension()
        N = lincomb.get_number_components()
        J = lincomb.get_number_packets()
        Cj = squeeze(lincomb.get_coefficients())
        eps = lincomb.get_eps()
        G = GradientHAWP()

        new_lincombs = [ LinearCombinationOfHAWPs(D, N, eps) for d in range(D) ]

        # Handle each wavepacket individually
        for j in range(J):
            packet = lincomb.get_wavepacket(j)
            grads = G.apply_gradient(packet, component=component)

            for d, grad in enumerate(grads):
                new_lincombs[d].add_wavepacket(grad, Cj[j])

        return new_lincombs
