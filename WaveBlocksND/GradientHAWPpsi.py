"""The WaveBlocks Project

Compute the action of the gradient operator applied to a new-kind Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, squeeze, dot, diag, real
from numpy.linalg import inv
from scipy import sqrt
from scipy.linalg import polar, eigh

from WaveBlocksND.WavepacketGradient import WavepacketGradient

__all__ = ["GradientHAWPpsi"]


class GradientHAWPpsi(WavepacketGradient):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x` applied to
    a new-kind Hagedorn wavepacket :math:`\Psi`.
    """

    def apply_gradient_component(self, wavepacket, component):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x`
        on the basis functions :math:`\psi(x)` of a component :math:`\Phi_i` of the
        new-kind Hagedorn wavepacket :math:`\Psi`.

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
        q, p, Q, P, _ = wavepacket.get_parameters(component=component)

        _, PA = polar(Q, side='left')
        EW, EV = eigh(real(PA))

        E = real(dot(P, inv(Q)))
        F1 = dot(E, dot(inv(EV.T), diag(EW)))
        F2 = 1.0j * dot(inv(EV.T), diag(1.0 / EW))
        Gb = F1 - F2
        Gf = F1 + F2

        coeffs = wavepacket.get_coefficients(component=component)

        # Prepare storage for new coefficients
        K = wavepacket.get_basis_shapes(component=component)
        Ke = K.extend()
        size = Ke.get_basis_size()
        cnew = zeros((size, D), dtype=complexfloating)

        # We implement the more efficient scatter type stencil here
        for k in K.get_node_iterator():
            # Central phi_i coefficient
            cnew[Ke[k], :] += squeeze(coeffs[K[k]] * p)

            # Backward neighbours phi_{i - e_d}
            nbw = Ke.get_neighbours(k, selection="backward")

            for d, nb in nbw:
                cnew[Ke[nb], :] += (sqrt(eps**2 / 2.0) *
                                    sqrt(k[d]) * coeffs[K[k]] *
                                    Gb[:, d])

            # Forward neighbours phi_{i + e_d}
            nfw = Ke.get_neighbours(k, selection="forward")

            for d, nb in nfw:
                cnew[Ke[nb], :] += (sqrt(eps**2 / 2.0) *
                                    sqrt(k[d] + 1.0) * coeffs[K[k]] *
                                    Gf[:, d])

        return (Ke, cnew)
