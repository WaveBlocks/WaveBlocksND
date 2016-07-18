"""The WaveBlocks Project

Compute the action of the position operator applied to a Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, conjugate, squeeze
from scipy import sqrt

from WaveBlocksND.WavepacketPosition import WavepacketPosition

__all__ = ["PositionHAWP"]


class PositionHAWP(WavepacketPosition):
    r"""This class implements the computation of the action of the
    position operator :math:`x` applied to a Hagedorn wavepacket :math:`\Psi`.
    """

    def __init__(self):
        r"""
        """
        pass


    def apply_position_component(self, wavepacket, component):
        r"""Compute the effect of the position operator :math:`x` on the basis functions :math:`\phi(x)`
        of a component :math:`\Phi_i` of the Hagedorn wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer.
        :return: Extended basis shape :math:`\mathfrak{\dot{K}}` and new coefficients :math:`c^\prime`
                 for component :math:`\Phi_i`. The coefficients are stored column-wise with
                 one column per dimension :math:`d`. The :math:`c^\prime` array is of shape
                 :math:`|\mathfrak{\dot{K}}| \times D`.
        """
        D = wavepacket.get_dimension()
        eps = wavepacket.get_eps()
        q, p, Q, P, S = wavepacket.get_parameters(component=component)
        Qbar = conjugate(Q)
        coeffs = wavepacket.get_coefficients(component=component)

        # Prepare storage for new coefficients
        K = wavepacket.get_basis_shapes(component=component)
        Ke = K.extend()
        size = Ke.get_basis_size()
        cnew = zeros((size, D), dtype=complexfloating)

        # We implement the more efficient scatter type stencil here
        for k in K.get_node_iterator():
            # Central phi_i coefficient
            cnew[Ke[k], :] += squeeze(coeffs[K[k]] * q)

            # Backward neighbours phi_{i - e_d}
            nbw = Ke.get_neighbours(k, selection="backward")

            for d, nb in nbw:
                cnew[Ke[nb], :] += sqrt(eps**2 / 2.0) * sqrt(k[d]) * coeffs[K[k]] * Qbar[:, d]

            # Forward neighbours phi_{i + e_d}
            nfw = Ke.get_neighbours(k, selection="forward")

            for d, nb in nfw:
                cnew[Ke[nb], :] += sqrt(eps**2 / 2.0) * sqrt(k[d] + 1.0) * coeffs[K[k]] * Q[:, d]

        return (Ke, cnew)
