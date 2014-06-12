"""The WaveBlocks Project

Compute some observables like norm, kinetic and potential energy
of Hagedorn wavepackets. This class implements the mixed case
where the bra does not equal the ket.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from functools import partial
from numpy import squeeze, sum

from Observables import Observables
from GradientHAWP import GradientHAWP

__all__ = ["ObservablesMixedHAWP"]


class ObservablesMixedHAWP(Observables):
    r"""This class implements observable computation for Hagedorn wavepackets :math:`\Psi`.
    This class implements the mixed case :math:`\langle \Psi | \cdot | \Psi^\prime \rangle`
    where the bra does not equal the ket.
    """

    # TODO: Support multi-component wavepackets

    def __init__(self, innerproduct=None):
        r"""Initialize a new :py:class:`ObservablesHAWP` instance for observable computation
        of Hagedorn wavepackets.

        :param innerproduct: An inner product for computing the integrals. The inner product is only used for
                             the computation of the potential energy :math:`\langle\Psi|V(x)|\Psi\rangle`
                             but not for the kinetic energy.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        # A innerproduct to compute the integrals
        if innerproduct is not None:
            self._innerproduct = innerproduct

        self._gradient = GradientHAWP()


    def set_innerproduct(self, innerproduct):
        r"""Set the innerproduct.

        :param innerproduct: An innerproduct for computing the integrals. The inner product is only used for
                             the computation of the potential energy :math:`\langle\Psi|V(x)|\Psi^\prime\rangle`
                             but not for the kinetic energy.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        self._innerproduct = innerproduct


    def overlap(self, pacbra, packet, component=None, summed=False):
        r"""Calculate the overlap :math:`\langle\Psi|\Psi^\prime\rangle` of the
        wavepackets :math:`\Psi` and :math:`\Psi^\prime`.

        :param pacbra: The wavepacket :math:`\Psi` of which we compute the overlap.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` of which we compute the overlap.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose norm is calculated.
                          The default value is ``None`` which means to compute the norms of all :math:`N` components.
        :type component: int or ``None``.
        :param summed: Whether to sum up the norms :math:`\langle\Phi_i|\Phi_i\rangle` of the
                       individual components :math:`\Phi_i`.
        :type summed: Boolean, default is ``False``.
        :type summed: Boolean, default is ``False``.
        :return: The norm of :math:`\Psi` or the norm of :math:`\Phi_i` or a list with the :math:`N`
                 norms of all components. Depending on the values of ``component`` and ``summed``.
        """
        return self._innerproduct.quadrature(pacbra, packet, component=0, summed=summed)


    def norm(self, wavepacket, component=None, summed=False):
        r"""
        """
        return self.overlap(wavepacket, wavepacket, component=component, summed=summed)


    def kinetic_overlap_energy(self, pacbra, packet, componentbra=None, componentket=None, summed=False):
        r"""Compute the kinetic energy :math:`E_{\text{kin}} := \langle\Psi|T|\Psi\rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param pacbra: The wavepacket :math:`\Psi` of which we compute the kinetic energy.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` of which we compute the kinetic energy.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose
                          kinetic energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param summed: Whether to sum up the kinetic energies :math:`E_i` of the individual
                       components :math:`\Phi_i`. Default is ``False``.
        :type summed: Boolean
        :return: A list with the kinetic energies of the individual components or the
                 overall kinetic energy of the wavepacket. (Depending on the optional arguments.)
        """
        if componentbra is None:
            Nbra = pacbra.get_number_components()
            componentsbra = xrange(Nbra)
        else:
            componentsbra = [componentbra]
        if componentket is None:
            Nket = packet.get_number_components()
            componentsket = xrange(Nket)
        else:
            componentsket = [componentket]

        ekin = []

        for nr in componentsbra:
            for nc in componentsket:
                #Kprimebra, cnewbra = self._gradient.apply_gradient(pacbra, nr, as_packet=True)
                #Kprimeket, cnewket = self._gradient.apply_gradient(packet, nc, as_packet=True)
                #ekin.append(0.5 * sum(sum(conjugate(cnewbra)*cnewket, axis=1), axis=0))
                gradpacbra = self._gradient.apply_gradient(pacbra, nr, as_packet=True)
                gradpacket = self._gradient.apply_gradient(packet, nc, as_packet=True)
                Q = self._innerproduct.quadrature(gradpacbra[0], gradpacket[0], summed=True)
                ekin.append(0.5*Q)

        if summed is True:
            ekin = sum(ekin)
        elif componentbra is not None or componentket is not None:
            # Do not return a list for specific single components
            ekin = ekin[0]

        return ekin


    def kinetic_energy(self, wavepacket, component=None, summed=False):
        r"""
        """
        return self.kinetic_overlap_energy(wavepacket, wavepacket, componentbra=component, componentket=component, summed=summed)


    def potential_overlap_energy(self, pacbra, packet, potential, component=None, summed=False):
        r"""Compute the potential energy :math:`E_{\text{pot}} := \langle\Psi|V|\Psi\rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param pacbra: The wavepacket :math:`\Psi` of which we compute the potential energy.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` of which we compute the potential energy.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param potential: The potential :math:`V(x)`. (Actually, not the potential object itself
                          but one of its ``V.evaluate_*`` methods.)
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose
                          potential energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param summed: Whether to sum up the potential energies :math:`E_i` of the individual
                       components :math:`\Phi_i`. Default is ``False``.
        :type summed: Boolean
        :return: A list with the potential energies of the individual components or the
                 overall potential energy of the wavepacket. (Depending on the optional arguments.)
        """
        Nbra = pacbra.get_number_components()
        Nket = packet.get_number_components()

        # TODO: Better take 'V' instead of 'V.evaluate_at' as argument?
        #f = partial(potential.evaluate_at, as_matrix=True)
        f = partial(potential, as_matrix=True)

        # Compute the brakets for each component
        if component is not None:
            Q = self._innerproduct.quadrature(pacbra, packet, operator=f, diag_component=component, eval_at_once=True)
        else:
            Q = self._innerproduct.quadrature(pacbra, packet, operator=f, eval_at_once=True)

        # And don't forget the summation in the matrix multiplication of 'operator' and 'ket'
        # TODO: Should this go inside the innerproduct?
        tmp = map(squeeze, Q)
        epot = [ sum(tmp[i*Nket:(i+1)*Nket]) for i in xrange(Nbra) ]

        if summed is True:
            epot = sum(epot)
        elif component is not None:
            # Do not return a list for specific single components
            epot = epot[0]

        return epot


    def potential_energy(self, wavepacket, potential, component=None, summed=False):
        r"""
        """
        return self.potential_overlap_energy(wavepacket, wavepacket, potential, component=component, summed=summed)
