"""The WaveBlocks Project

Compute some observables like kinetic and potential energy
of a Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, vstack, conjugate, dot, squeeze, sum
from scipy import sqrt
from scipy.linalg import norm

from Observables import Observables

__all__ = ["ObservablesHAWP"]


class ObservablesHAWP(Observables):
    r"""This class implements observable computation for Hagedorn wavepackets.
    """

    def __init__(self, quadrature=None):
        r"""
        """
        # A quadrature to compute the integrals
        if quadrature is not None:
            self._quadrature = quadrature


    def set_quadrature(self, quadrature):
        r"""
        """
        self._quadrature = quadrature


    def apply_gradient(self, wavepacket, component):
        r"""Compute the effect of the operator :math:`-i \varepsilon^2 \nabla_x` on the basis
        functions :math:`\phi` of a component :math:`\Phi_i` of the Hagedorn wavepacket :math:`\Psi`.

        :param wavepacket:
        :param component:
        :return: Extended basis shape :math:`\mathcal{\dot{K}}` and new coefficients :math:`c^\prime`.
        """
        # TODO: Maybe move this method into the HAWP class?
        D = wavepacket.get_dimension()
        K = wavepacket.get_basis_shape(component=component)

        eps = wavepacket.get_eps()
        q, p, Q, P, S = wavepacket.get_parameters(component=component)
        Pbar = conjugate(P)
        coeffs = wavepacket.get_coefficients(component=component)

        size = K.get_basis_size(extended=True)
        cnew = zeros((size,D), dtype=complexfloating)

        # We implement the less efficient gather type stencil here
        I = K.get_node_iterator(extended=True)

        for i in I:
            # Current phi_i coefficient
            if i in K:
                ccur = coeffs[K[i]]
            else:
                ccur = 0.0

            # Backward neighbours phi_i - <e>
            cbw = zeros((D,1), dtype=complexfloating)
            nbw = K.get_neighbours(i, selection="backward", extended=False)

            for d, nb in nbw:
                j = K[nb]
                cbw[d] = coeffs[j] * sqrt(i[d])

            # Forward neighbours phi_i + <e>
            nfw = K.get_neighbours(i, selection="forward", extended=False)
            cfw = zeros((D,1), dtype=complexfloating)

            for d, nb in nfw:
                j = K[nb]
                cfw[d] = coeffs[j] * sqrt(i[d] + 1.0)

            # Compute parts and assemble
            res = sqrt(eps**2/2.0) * (dot(P, cfw) + dot(Pbar, cbw)) + p * ccur
            cnew[K[i],:] = squeeze(res)

        return (K, cnew)


    def kinetic_energy(self, wavepacket, component=None, summed=False):
        r"""Calculate the kinetic energy :math:`E_{\text{kin}} := \langle\Psi|T|\Psi\rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` of which we compute the kinetic energy.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component:
        :param summed: Whether to sum up the kinetic energies :math:`E_i` of the individual
                       components :math:`\Phi_i`. Default is `False`.
        :return: A list with the kinetic energies of the individual components or the
                 overall kinetic energy of the wavepacket. (Depending on the optional arguments.)
        """
        N = wavepacket.get_number_components()

        ekin = []

        for n in xrange(N):
            Kprime, cnew = self.apply_gradient(wavepacket, n)
            ekin.append(0.5 * sum(sum(conjugate(cnew)*cnew, axis=1), axis=0))

        if summed is True:
            ekin = sum(ekin)

        return ekin


    def potential_energy(self, wavepacket, potential, summed=False):
        r"""Calculate the potential energy :math:`E_{\text{pot}} := \langle\Psi|V|\Psi\rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` of which we compute the potential energy.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param potential: The potential :math:`V(x)`.
        :param summed: Whether to sum up the potential energies :math:`E_i` of the individual
                       components :math:`\Phi_i`. Default is `False`.
        :return: A list with the potential energies of the individual components or the
                 overall potential energy of the wavepacket. (Depending on the optional arguments.)
        """
        pass
