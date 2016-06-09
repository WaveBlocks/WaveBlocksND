"""The WaveBlocks Project

Compute the action of the gradient operator applied to a new Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, conjugate, squeeze, dot, diag, real
from numpy.linalg import inv
from scipy import sqrt
from scipy.linalg import polar, eigh

from WaveBlocksND.WavepacketGradient import WavepacketGradient

__all__ = ["GradientHAWPnew"]


class GradientHAWPnew(WavepacketGradient):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x` applied to
    a Hagedorn wavepacket :math:`\Psi`.
    """

    def __init__(self):
        r"""
        """
        pass


    def apply_gradient_component(self, wavepacket, component):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x` on the basis
        functions :math:`\phi(x)` of a component :math:`\Phi_i` of the Hagedorn wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer.
        :return: Extended basis shape :math:`\mathfrak{\dot{K}}` and new coefficients :math:`c^\prime`
                 for component :math:`\Phi_i`. The coefficients are stored column-wise with
                 one column per dimension :math:`d`. The :math:`c^\prime` array is of shape
                 :math:`|\mathfrak{\dot{K}}| \times D`.
        """
        # TODO: Remove
        if not wavepacket._new:
            raise ValueError("Old style wavepacket in new style gradient operator!")

        D = wavepacket.get_dimension()
        eps = wavepacket.get_eps()
        q, p, Q, P, S = wavepacket.get_parameters(component=component)

        _, PA = polar(Q, side='left')
        EW, EV = eigh(real(PA))

        C = -1.0j * conjugate(conjugate(dot(P, inv(Q))) + dot(P, inv(Q))).T
        Gb = 0.5 * dot(C, dot(inv(EV.T), diag(EW))) - dot(inv(EV.T), diag(1.0 / EW))
        Gf = 0.5 * dot(C, dot(inv(EV.T), diag(EW))) + dot(inv(EV.T), diag(1.0 / EW))

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
                cnew[Ke[nb], :] += (1.0j * sqrt(eps**2 / 2.0) *
                                    sqrt(k[d]) * coeffs[K[k]] *
                                    Gb[:, d])

            # Forward neighbours phi_{i + e_d}
            nfw = Ke.get_neighbours(k, selection="forward")

            for d, nb in nfw:
                cnew[Ke[nb], :] += (1.0j * sqrt(eps**2 / 2.0) *
                                    sqrt(k[d] + 1.0) * coeffs[K[k]] *
                                    Gf[:, d])

        return (Ke, cnew)
