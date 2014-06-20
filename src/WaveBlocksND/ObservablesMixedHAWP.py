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
    r"""This class implements the mixed case observable computation
    :math:`\langle \Psi | \cdot | \Psi^\prime \rangle` for Hagedorn
    wavepackets :math:`\Psi` where the bra :math:`\Psi` does not equal
    the ket :math:`\Psi^\prime`.

    .. warning:: Wavepackets consisting of multiple components are not yet supported.
    """

    # TODO: Support multi-component wavepackets

    def __init__(self, innerproduct=None):
        r"""Initialize a new :py:class:`ObservablesMixedHAWP` instance for observable computation of
        Hagedorn wavepackets.

        :param innerproduct: An inner product for computing the integrals. The inner product is used
                             for the computation of all brakets
                             :math:`\langle \Psi | \cdot | \Psi^\prime \rangle`.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        # A innerproduct to compute the integrals
        if innerproduct is not None:
            self._innerproduct = innerproduct

        self._gradient = GradientHAWP()


    def set_innerproduct(self, innerproduct):
        r"""Set the innerproduct.

        :param innerproduct: An inner product for computing the integrals. The inner product is used
                             for the computation of all brakets
                             :math:`\langle \Psi | \cdot | \Psi^\prime \rangle`.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        self._innerproduct = innerproduct


    def overlap(self, pacbra, packet, componentbra=None, componentket=None, summed=False):
        r"""Calculate the overlap :math:`\langle \Psi | \Psi^\prime \rangle` of the wavepackets
        :math:`\Psi` and :math:`\Psi^\prime`.

        :param pacbra: The wavepacket :math:`\Psi` which takes part in the overlap integral.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` which takes part in the overlap integral.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param componentbra: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`
                             whose overlap is computed. The default value is ``None`` which means
                             to compute the overlaps with all :math:`N` components of :math:`\Psi`.
        :type componentbra: Integer or ``None``.
        :param componentket: The index :math:`i` of the component :math:`\Phi_i^\prime` of :math:`\Psi^\prime`
                             whose overlap is computed. The default value is ``None`` which means
                             to compute the overlaps with all :math:`N^\prime` components of :math:`\Psi^\prime`.
        :type componentket: Integer or ``None``.
        :param summed: Whether to sum up the overlaps :math:`\langle \Phi_i | \Phi_i^\prime \rangle`
                       of the individual components :math:`\Phi_i` and :math:`\Phi_i^\prime`.
        :type summed: Boolean, default is ``False``.
        :return: The overlap of :math:`\Psi` with :math:`\Psi^\prime` or the overlap of :math:`\Phi_i`
                 with :math:`\Phi_i^\prime` or a list with the :math:`N` overlaps of all components.
                 (Depending on the optional arguments.)
        """
        return self._innerproduct.quadrature(pacbra, packet, component=0, summed=summed)


    def norm(self, wavepacket, component=None, summed=False):
        r"""Calculate the :math:`L^2` norm :math:`\langle \Psi | \Psi \rangle` of the wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` of which we compute the norm.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose norm is computed.
                          The default value is ``None`` which means to compute the norms of all :math:`N` components.
        :type component: int or ``None``.
        :param summed: Whether to sum up the norms :math:`\langle \Phi_i | \Phi_i \rangle` of the
                       individual components :math:`\Phi_i`.
        :type summed: Boolean, default is ``False``.
        :return: The norm of :math:`\Psi` or the norm of :math:`\Phi_i` or a list with the :math:`N`
                 norms of all components. (Depending on the optional arguments.)

        .. note:: This method just expands to a call of the :py:meth:`ObservablesMixedHAWP.overlap`
                  method. Better use :py:meth:`ObservablesHAWP.norm`.
        """
        return self.overlap(wavepacket, wavepacket, componentbra=component, componentket=component, summed=summed)


    def kinetic_overlap_energy(self, pacbra, packet, componentbra=None, componentket=None, summed=False):
        r"""Compute the kinetic energy overlap :math:`\langle \Psi | T | \Psi^\prime \rangle`
        of the different components :math:`\Phi_i` and :math:`\Phi_i^\prime` of the
        wavepackets :math:`\Psi` and :math:`\Psi^\prime`.

        :param pacbra: The wavepacket :math:`\Psi` which takes part in the kinetic energy integral.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` which takes part in the kinetic energy integral.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param componentbra: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`
                             which takes part in the kinetic energy integral. If set to ``None`` the
                             computation is performed for all :math:`N` components of :math:`\Psi`.
        :type componentbra: Integer or ``None``.
        :param componentket: The index :math:`i` of the component :math:`\Phi_i^\prime` of :math:`\Psi^\prime`
                             which takes part in the kinetic energy integral. If set to ``None`` the
                             computation is performed for all :math:`N^\prime` components of :math:`\Psi^\prime`.
        :type componentket: Integer or ``None``.
        :param summed: Whether to sum up the kinetic energies :math:`E_i` of the individual
                       components :math:`\Phi_i` and :math:`\Phi_i^\prime`.
        :type summed: Boolean, default is ``False``.
        :return: A list of the kinetic energy overlap integrals of the individual components or
                 the overall kinetic energy overlap of the wavepackets. (Depending on the optional arguments.)
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
                gradpacbra = self._gradient.apply_gradient(pacbra, nr, as_packet=True)
                gradpacket = self._gradient.apply_gradient(packet, nc, as_packet=True)
                Q = [self._innerproduct.quadrature(gpb, gpk, summed=True) for gpb, gpk in zip(gradpacbra, gradpacket)]
                ekin.append(0.5*sum(Q))

        if summed is True:
            ekin = sum(ekin)
        elif componentbra is not None or componentket is not None:
            # Do not return a list for specific single components
            ekin = ekin[0]

        return ekin


    def kinetic_energy(self, wavepacket, component=None, summed=False):
        r"""Compute the kinetic energy :math:`E_{\text{kin}} := \langle \Psi | T | \Psi \rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` of which we compute the kinetic energy.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose
                          kinetic energy we compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param summed: Whether to sum up the kinetic energies :math:`E_i` of the individual
                       components :math:`\Phi_i`.
        :type summed: Boolean, default is ``False``.
        :return: A list of the kinetic energies of the individual components or the
                 overall kinetic energy of the wavepacket. (Depending on the optional arguments.)

        .. note:: This method just expands to a call of the :py:meth:`ObservablesMixedHAWP.kinetic_overlap_energy`
                  method. Better use :py:meth:`ObservablesHAWP.kinetic_energy`.
        """
        return self.kinetic_overlap_energy(wavepacket, wavepacket, componentbra=component, componentket=component, summed=summed)


    def potential_overlap_energy(self, pacbra, packet, potential, componentbra=None, componentket=None, summed=False):
        r"""Compute the potential energy overlap :math:`\langle \Psi | V(x) | \Psi^\prime \rangle`
        of the different components :math:`\Phi_i` and :math:`\Phi_i^\prime` of the
        wavepackets :math:`\Psi` and :math:`\Psi^\prime`.

        :param pacbra: The wavepacket :math:`\Psi` which takes part in the potential energy integral.
        :type pacbra: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param packet: The wavepacket :math:`\Psi^\prime` which takes part in the potential energy integral.
        :type packet: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param potential: The potential :math:`V(x)`. (Actually, not the potential object itself
                          but one of its ``V.evaluate_*`` methods.)
        :param componentbra: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`
                             which takes part in the potential energy integral. If set to ``None`` the
                             computation is performed for all :math:`N` components of :math:`\Psi`.
        :type componentbra: Integer or ``None``.
        :param componentket: The index :math:`i` of the component :math:`\Phi_i^\prime` of :math:`\Psi^\prime`
                             which takes part in the potential energy integral. If set to ``None`` the
                             computation is performed for all :math:`N^\prime` components of :math:`\Psi^\prime`.
        :type componentket: Integer or ``None``.
        :param summed: Whether to sum up the potential energies :math:`E_i` of the individual
                       components :math:`\Phi_i` and :math:`\Phi_i^\prime`.
        :type summed: Boolean, default is ``False``.
        :return: A list of the potential energy overlap integrals of the individual components or
                 the overall potential energy overlap of the wavepackets. (Depending on the optional arguments.)
        """
        Nbra = pacbra.get_number_components()
        Nket = packet.get_number_components()

        if componentbra is None:
            componentsbra = xrange(Nbra)
        else:
            componentsbra = [componentbra]
        if componentket is None:
            componentsket = xrange(Nket)
        else:
            componentsket = [componentket]

        # TODO: Fix
        component = None

        # TODO: Better take 'V' instead of 'V.evaluate_at' as argument?
        #f = partial(potential.evaluate_at, as_matrix=True)
        f = partial(potential, as_matrix=True)

        # Compute the brakets for each component
        if component is not None:
            Q = self._innerproduct.quadrature(pacbra, packet, operator=f, diag_component=component, eval_at_once=True)
            Q = [squeeze(Q)]
        else:
            Q = self._innerproduct.quadrature(pacbra, packet, operator=f, eval_at_once=True)
            Q = map(squeeze, Q)

        # And don't forget the summation in the matrix multiplication of 'operator' and 'ket'
        # TODO: Should this go inside the innerproduct?
        epot = [ sum(Q[i*Nket:(i+1)*Nket]) for i in xrange(Nbra) ]

        if summed is True:
            epot = sum(epot)
        elif component is not None:
            # Do not return a list for specific single components
            epot = epot[0]

        return epot


    def potential_energy(self, wavepacket, potential, component=None, summed=False):
        r"""Compute the potential energy :math:`E_{\text{pot}} := \langle \Psi | V(x) | \Psi \rangle`
        of the different components :math:`\Phi_i` of the wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` of which we compute the potential energy.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param potential: The potential :math:`V(x)`. (Actually, not the potential object itself
                          but one of its ``V.evaluate_*`` methods.)
        :param component: The index :math:`i` of the component :math:`\Phi_i` whose
                          potential energy we compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param summed: Whether to sum up the potential energies :math:`E_i` of the individual
                       components :math:`\Phi_i`.
        :type summed: Boolean, default is ``False``.
        :return: A list of the potential energies of the individual components or the
                 overall potential energy of the wavepacket. (Depending on the optional arguments.)

        .. note:: This method just expands to a call of the :py:meth:`ObservablesMixedHAWP.potential_overlap_energy`
                  method. Better use :py:meth:`ObservablesHAWP.potential_energy`.
        """
        return self.potential_overlap_energy(wavepacket, wavepacket, potential, componentbra=component, componentket=component, summed=summed)
