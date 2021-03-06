"""The WaveBlocks Project

This file contains the Magnus propagator class for homogeneous wavepackets.

@author: V. Gradinaru
@copyright: Copyright (C) 2012, 2014, 2015, 2016 V. Gradinaru, R. Bourquin
@license: Modified BSD License
"""

from functools import partial
from numpy import dot, eye, atleast_2d, sqrt
from numpy.linalg import inv, det

from WaveBlocksND.Propagator import Propagator
from WaveBlocksND.BlockFactory import BlockFactory
from WaveBlocksND.SplittingParameters import SplittingParameters
from WaveBlocksND.ComplexMath import cont_angle

__all__ = ["MagnusPropagator"]


class MagnusPropagator(Propagator, SplittingParameters):
    r"""This class can numerically propagate given initial values :math:`\Psi` in
    a potential :math:`V(x)`. The propagation is done for a given set of homogeneous
    Hagedorn wavepackets neglecting interaction."""

    def __init__(self, parameters, potential, packets=[]):
        r"""Initialize a new :py:class:`MagnusPropagator` instance.

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

        self._a, self._b = self.build(parameters["splitting_method"])


    def __str__(self):
        r"""Prepare a printable string representing the :py:class:`MagnusPropagator` instance."""
        return "Homogeneous Magnus propagator for " + str(self._number_components) + " components.\n"


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


    # TODO: Consider removing this, duplicate
    def get_number_components(self):
        r""":return: The number :math:`N` of components :math:`\Phi_i` of :math:`\Psi`.
        """
        return self._number_components


    def get_wavepackets(self, packet=None):
        r"""Return the wavepackets :math:`\{\Psi_i\}_i` that take part in the time propagation by the
        current :py:class:`MagnusPropagator` instance.

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


    def propagate(self):
        r"""Given a wavepacket :math:`\Psi` at time :math:`t` compute the propagated
        wavepacket at time :math:`t + \tau`. We perform exactly one timestep of size
        :math:`\tau` here. This propagation is done for all packets in the list
        :math:`\{\Psi_i\}_i` and neglects any interaction between two packets.
        The semiclassical propagation scheme is used.

        The Gauss-Legendre rule on the interval :math:`[-1,1]` has two nodes
        :math:`n_1 = -\frac{1}{\sqrt{3}}` and
        :math:`n_2 = \frac{1}{\sqrt{3}}`.
        The integration time step is :math:`[0, 1] dt`, hence we rescale the nodes
        onto the interval :math:`[0, 1]` to get
        :math:`n_1^{\prime} = \frac{1}{2} - \frac{1}{2\sqrt{3}} =  \frac{1}{2} - \frac{\sqrt{3}}{6}` and
        :math:`n_2^{\prime} = \frac{1}{2} + \frac{1}{2\sqrt{3}} =  \frac{1}{2} + \frac{\sqrt{3}}{6}`.
        The two time points :math:`h_1` and :math:`h_2` are then given by
        :math:`h_1 = n_1^{\prime} = \frac{1}{2} - \frac{1}{2\sqrt{3}}` and
        :math:`h_2 = n_2^{\prime} - n_1^{\prime} = \frac{1}{\sqrt{3}}`.
        Given a (matrix) differential equation:

        .. math::

            y^{\prime}(t) = A(t) y(t)

        with :math:`t \geq 0`. In each timestep of size :math:`h = dt` we compute:

        .. math::

            A_1 &= A \left( t_n + \left(\frac{1}{2} - \frac{\sqrt{3}}{6} \right) h \right) \\
            A_2 &= A \left( t_n + \left(\frac{1}{2} + \frac{\sqrt{3}}{6} \right) h \right)

        from which we get:

        .. math::

            \sigma_n = \frac{1}{2} h (A_1 + A_2) + \frac{\sqrt{3}}{12} h^2 [A_2, A_1]

        and then we find for the solution:

        .. math::

            y_{n+1} = e^{\sigma_n} y_n

        More details can be found in [#]_ and [#]_ and especially formula 2.9 in [#]_.

        .. [#] S. Blanes and P.C. Moan, "Fourth- and sixth-order commutator-free Magnus integrators for linear and non-linear dynamical systems",
               Applied Numerical Mathematics, volume 56 number 12 (2006) 1519-1537.

        .. [#] S. Blanes and F. Casas and J. Ros, "Improved high order integrators based on the Magnus expansion",
               BIT Numerical Mathematics, volume 40 (1999) 434-450.

        .. [#] A. Iserles, A. Marthinsen and S.P. Norsett, "On the Implementation of the Method of Magnus Series for Linear Differential Equations",
               BIT Numerical Mathematics, volume 39 number 2 (1999) 281-304.
        """
        # Cache some parameter values
        dt = self._dt
        a = self._a
        b = self._b

        # Propagate all packets
        for packet, leading_chi in self._packets:
            eps = packet.get_eps()

            # Propagate until c1 * dt
            h1 = (0.5 - sqrt(3.0) / 6.0) * dt

            # Inner time step
            nrinnersteps = self._parameters.get("innersteps", h1**0.5 * eps**(-3.0 / 8.0))
            nrlocalsteps1 = max(1, 1 + int(nrinnersteps))
            self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, h1], nrlocalsteps1, [packet], [packet, leading_chi])

            # Build a first matrix here with the current parameters of the wavepacket
            innerproduct = packet.get_innerproduct()
            A1 = -1.0j / eps**2 * innerproduct.build_matrix(packet, operator=partial(self._potential.evaluate_local_remainder_at, diagonal_component=leading_chi))

            # Propagate until c2 * dt
            h2 = 1.0 / sqrt(3.0) * dt

            # Inner time step
            nrinnersteps = self._parameters.get("innersteps", h2**0.5 * eps**(-3.0 / 8.0))
            nrlocalsteps2 = max(1, 1 + int(nrinnersteps))
            self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, h2], nrlocalsteps2, [packet], [packet, leading_chi])

            # Build a second matrix here with the current parameters of the wavepacket
            innerproduct = packet.get_innerproduct()
            A2 = -1.0j / eps**2 * innerproduct.build_matrix(packet, operator=partial(self._potential.evaluate_local_remainder_at, diagonal_component=leading_chi))

            # Combine A1 and A2 and build the F matrix for Magnus of 4-th order split
            F = 0.5 * dt * (A1 + A2) + sqrt(3.0) / 12.0 * dt**2 * (dot(A2, A1) - dot(A1, A2))

            # Propagate the coefficients
            coefficients = packet.get_coefficient_vector()
            coefficients = self._matrix_exponential(F, coefficients, 1.0)
            packet.set_coefficient_vector(coefficients)

            # Propagate until dt to finish the current timestep
            self.intsplit(self._propkin, self._proppotquad, a, b, [0.0, h1], nrlocalsteps1, [packet], [packet, leading_chi])
