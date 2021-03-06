"""The WaveBlocks Project

This file contains a propagator class for homogeneous wavepackets.

@author: V. Gradinaru
@copyright: Copyright (C) 2013, 2014, 2015 V. Gradinaru, R. Bourquin
@license: Modified BSD License
"""

from functools import partial
from numpy import dot, eye, atleast_2d
from numpy.linalg import inv, det

from WaveBlocksND.Propagator import Propagator
from WaveBlocksND.BlockFactory import BlockFactory
from WaveBlocksND.SplittingParameters import SplittingParameters
from WaveBlocksND.ProcessingSplittingParameters import ProcessingSplittingParameters
from WaveBlocksND.ComplexMath import cont_angle

__all__ = ["Pre764scPropagator"]


class Pre764scPropagator(Propagator, SplittingParameters):
    r"""This class can numerically propagate given initial values :math:`\Psi` in
    a potential :math:`V(x)`. The propagation is done for a given set of homogeneous
    Hagedorn wavepackets neglecting interaction. It uses the preprocessor 764 method,
    see Blanes, Casas, Ros 2000, table IV and the idea of the semiclassical splitting
    """

    def __init__(self, parameters, potential, packets=[]):
        r"""Initialize a new :py:class:`SemiclassicalPropagator` instance.

        :param parameters: A :py:class:`ParameterProvider` instance containing at least
                           the key ``dt`` for providing the timestep :math:`\tau`.
        :type parameters: A :py:class:`ParameterProvider` instance
        :param potential: The potential :math:`V(x)` the wavepacket :math:`\Psi` feels during the time propagation.
        :param packet: The initial homogeneous Hagedorn wavepacket :math:`\Psi` we propagate in time.

        :raises ValueError: If the number of components of :math:`\Psi` does not match
                            the number of energy levels :math:`\lambda_i` of the potential.
        """
        # The potential :math:`V(x)` the packet(s) feel.
        self._potential = potential

        # Number :math:`N` of components the wavepacket :math:`\Psi` has got.
        self._number_components = self._potential.get_number_components()
        self._dimension = self._potential.get_dimension()

        # A list of Hagedorn wavepackets :math:`\Psi` together with some codata
        # like the leading component :math:`\chi` which is the index of the eigenvalue
        # :math:`\lambda_\chi` of the potential :math:`V` that is responsible for
        # propagating the Hagedorn parameters.
        # TODO: We assume a list of (packet, leading_component) tuples here. Generalize tuples to dicts!
        # TODO: Do not use a list but better use a hashtable by packet IDs?
        self._packets = packets[:]

        # Keep a reference to the parameter provider instance
        self._parameters = parameters

        self._dt = self._parameters["dt"]

        # The relative mass scaling matrix M
        if "mass_scaling" in self._parameters:
            self._M = atleast_2d(self._parameters["mass_scaling"])
            assert self._M.shape == (self._dimension, self._dimension)
            self._Minv = inv(self._M)
        else:
            # No mass matrix given. Scale all masses equally
            self._M = eye(self._dimension)
            self._Minv = self._M

        # Decide about the matrix exponential algorithm to use
        self.__dict__["_matrix_exponential"] = BlockFactory().create_matrixexponential(parameters)

        # Precalculate the potential splittings needed
        self._prepare_potential()
        self._a, self._b = SplittingParameters().build(parameters["splitting_method"])
        self._innerorder = SplittingParameters().order(parameters["splitting_method"])
        self._A, self._B, self._Y, self._Z = ProcessingSplittingParameters().build("BCR764")


    def __str__(self):
        r"""Prepare a printable string representing the :py:class:`SemiclassicalPropagator` instance."""
        return "Homogeneous Semiclassical propagator for " + str(self._number_components) + " components.\n"


    def _prepare_potential(self):
        """Precalculate the potential splittings needed
        """
        for chi in set([p[1] for p in self._packets]):
            self._potential.calculate_local_quadratic(diagonal_component=chi)
            self._potential.calculate_local_remainder(diagonal_component=chi)


    def add_wavepacket(self, packet):
        r"""Add a new wavepacket :math:`\Psi` to the list of propagated wavepackets.

        :param packet: The new wavepacket :math:`\Psi` and its leading component :math:`\chi \in [0,N-1]`.
        :type packet: A tuple :math:`(\Psi, \chi)` with :math:`\Psi` a :py:class:`HagedornWavepacket`
                      instance and :math:`\chi` an integer.
        """
        self._packets.append(tuple(packet))
        self._prepare_potential()


    # TODO: Consider removig this, duplicate
    def get_number_components(self):
        r""":return: The number :math:`N` of components :math:`\Phi_i` of :math:`\Psi`.
        """
        return self._number_components


    def get_wavepackets(self, packet=None):
        r"""Return the wavepackets :math:`\{\Psi_i\}_i` that take part in the time propagation by the
        current :py:class:`SemiclassicalPropagator` instance.

        :param packet: The index :math:`i` (in this list) of a single packet :math:`\Psi_i` that is
                       to be returned. If set to ``None`` (default) return the full list with all packets.
        :type packet: Integer or ``None``
        :return: A list of :py:class:`HagedornWavepacket` instances or a single instance.
        """
        # TODO: Does not return leading components. Add this if needed somewhere.
        if packet is None:
            return [p[0] for p in self._packets]
        else:
            return self._packets[packet][0]


    def set_wavepackets(self, packetlist):
        """Set the list :math:`\{\Psi_i\}_i` of wavepackets that the propagator will propagate.

        :param packetlist: A list of new wavepackets :math:`\Psi_i` and their
                           leading components :math:`\chi_i` to propagate.
        :type packetlist: A list of :math:`(\Psi_i, \chi_i)` tuples.
        """
        self._packets = packetlist[:]


    def _propkin(self, h, packet):
        """Do a kinetic step of size h.
        """
        Mi = self._Minv
        key = ("q", "p", "Q", "P", "S", "adQ")
        q, p, Q, P, S, adQ = packet.get_parameters(key=key)
        q = q + h * dot(Mi, p)
        Q = Q + h * dot(Mi, P)
        S = S + 0.5 * h * dot(p.T, dot(Mi, p))
        adQn = cont_angle(det(Q), reference=adQ)[0]
        packet.set_parameters((q, p, Q, P, S, adQn), key=key)


    def _proppotquad(self, h, packet, leading_chi):
        """Do a potential step of size h with the local quadratic part.
        """
        q, p, Q, P, S = packet.get_parameters()
        V = self._potential.evaluate_local_quadratic_at(q, diagonal_component=leading_chi)
        p = p - h * V[1]
        P = P - h * dot(V[2], Q)
        S = S - h * V[0]
        packet.set_parameters((q, p, Q, P, S))


    def pre_propagate(self):
        r"""Given the wavefunction :math:`\psi` at initial time :math:`t_0`,
        perform some computations exactly once before running the ordinary
        time propagation.
        """
        # Cache some parameter values
        dt = self._dt
        a = self._a
        b = self._b
        # Preprocessor
        Z = self._Z
        Y = self._Y
        v = Z.shape[0]
        # Apply preprocessor step
        for packet, leading_chi in self._packets:
            eps = packet.get_eps()

            # Inner time step (fit to third term: eps^3 dt^4)
            r = self._innerorder
            alpha = 3.0
            beta = 4.0
            defaultinnerstep = (dt**(r - beta) / eps**(alpha + 2.0))**(1.0 / r)
            nrinnersteps = self._parameters.get("innersteps", defaultinnerstep)
            nrlocalsteps = max(1, 1 + int(nrinnersteps))

            # Splitting
            for j in range(v):
                # Step with Abig
                self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, -Z[j] * dt], nrlocalsteps, [packet], [packet, leading_chi])

                # Step with Beps
                # Do a potential step with the local non-quadratic taylor remainder
                innerproduct = packet.get_innerproduct()
                F = innerproduct.build_matrix(packet, operator=partial(self._potential.evaluate_local_remainder_at, diagonal_component=leading_chi))
                coefficients = packet.get_coefficient_vector()
                coefficients = self._matrix_exponential(F, coefficients, 1.0j * Y[j] * dt / eps**2)
                packet.set_coefficient_vector(coefficients)


    def post_propagate(self):
        r"""Given the wavefunction :math:`\psi` at final time :math:`T`,
        perform some computations exactly once after running the ordinary
        time propagation.
        """
        # Cache some parameter values
        dt = self._dt
        a = self._a
        b = self._b
        # Postprocessor
        Z = self._Z
        Y = self._Y
        v = Z.shape[0]
        # Apply postprocessor step
        for packet, leading_chi in self._packets:
            eps = packet.get_eps()

            # Inner time step (fit to third term: eps^3 dt^4)
            r = self._innerorder
            alpha = 3.0
            beta = 4.0
            defaultinnerstep = (dt**(r - beta) / eps**(alpha + 2.0))**(1.0 / r)
            nrinnersteps = self._parameters.get("innersteps", defaultinnerstep)
            nrlocalsteps = max(1, 1 + int(nrinnersteps))

            # Splitting
            for j in reversed(range(v)):
                # Step with Beps
                # Do a potential step with the local non-quadratic taylor remainder
                innerproduct = packet.get_innerproduct()
                F = innerproduct.build_matrix(packet, operator=partial(self._potential.evaluate_local_remainder_at, diagonal_component=leading_chi))
                coefficients = packet.get_coefficient_vector()
                coefficients = self._matrix_exponential(F, coefficients, -1.0j * Y[j] * dt / eps**2)
                packet.set_coefficient_vector(coefficients)

                # Step with Abig
                self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, Z[j] * dt], nrlocalsteps, [packet], [packet, leading_chi])


    def propagate(self):
        r"""Given a wavepacket :math:`\Psi` at time :math:`t` compute the propagated
        wavepacket at time :math:`t + \tau`. We perform exactly one timestep of size
        :math:`\tau` here. This propagation is done for all packets in the list
        :math:`\{\Psi_i\}_i` and neglects any interaction between two packets.

        More details can be found in [#]_.

        .. [#] S. Blanes, R. Bourquin and V. Gradinaru, "Raising the Order of Convergence in the Semiclassical Splitting".
        """
        # Cache some parameter values
        dt = self._dt
        a = self._a
        b = self._b
        # Kernel Pattern ABA
        Beps = self._A
        Abig = self._B
        u = Beps.shape[0]

        # Propagate all packets
        for packet, leading_chi in self._packets:
            eps = packet.get_eps()

            # Inner time step (fit to third term: eps^3 dt^4)
            r = self._innerorder
            alpha = 3.0
            beta = 4.0
            defaultinnerstep = (dt**(r - beta) / eps**(alpha + 2.0))**(1.0 / r)
            nrinnersteps = self._parameters.get("innersteps", defaultinnerstep)
            nrlocalsteps = max(1, 1 + int(nrinnersteps))

            # Splitting
            for j in range(u):
                # Step with Beps
                # Avoid expensive computation if coefficient is zero
                if Beps[j] != 0.0:
                    innerproduct = packet.get_innerproduct()
                    F = innerproduct.build_matrix(packet, operator=partial(self._potential.evaluate_local_remainder_at, diagonal_component=leading_chi))
                    coefficients = packet.get_coefficient_vector()
                    coefficients = self._matrix_exponential(F, coefficients, -1.0j * Beps[j] * dt / eps**2)
                    packet.set_coefficient_vector(coefficients)

                # Step with Abig
                self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, Abig[j] * dt], nrlocalsteps, [packet], [packet, leading_chi])
