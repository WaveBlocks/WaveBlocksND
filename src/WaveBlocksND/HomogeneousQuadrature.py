"""The WaveBlocks Project

This file contains code for the homogeneous quadrature of Hagedorn wavepackets.
The class defined here can compute brakets, inner products and expectation
values and compute the matrix elements of an arbitrary operator.

@author: R. Bourquin
@copyright: Copyright (C) 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, ones, complexfloating, sum, cumsum, squeeze, conjugate, dot, outer, transpose
from scipy.linalg import sqrtm, inv

from Quadrature import Quadrature


class HomogeneousQuadrature(Quadrature):

    def __init__(self, QR=None):
        if QR is not None:
            self.set_qr(QR)
        else:
            self._QR = None


    def __str__(self):
        return "Homogeneous quadrature using a " + str(self.QR)


    def transform_nodes(self, Pi, eps, QR=None):
        r"""Transform the quadrature nodes :math:`\gamma` such that they fit the given wavepacket :math:`\Phi\left[\Pi\right]`.

        :param Pi: The parameter set :math:`\Pi` of the wavepacket.
        :param eps: The value of :math:`\varepsilon` of the wavepacket.
        :param QR: An optional quadrature rule :math:`(\gamma, \omega)` providing the nodes. If not given
                   the internal quadrature rule will be used.
        """
        if QR is None:
            QR = self._QR

        q, p, Q, P, S = Pi

        # TODO: Try to get rid of inverses by using LU
        Q0 = inv(dot(Q, conjugate(transpose(Q))))
        Qs = inv(sqrtm(Q0))

        # TODO: Avoid sqrtm computation, use svd
        #U, S, Vh = svd(...)
        #Sinv = diagsvd(s)
        #Qs = V Sinv Vh

        nodes = q + eps * dot(Qs, QR.get_nodes())
        return nodes.copy()


    def quadrature(self, packet, operator=None, summed=False, component=None, diag_component=None):
        r"""Performs the quadrature of :math:`\langle\Psi|f|\Psi\rangle` for a general
        function :math:`f(x)` with :math:`x \in \mathbb{R}^D`.


        :param packet: The wavepacket :math:`\Psi`.
        :param operator: A real-valued function :math:`f(x): \mathbb{R}^D \rightarrow \mathbb{R}^{N \times N}`.
        :param summed: Whether to sum up the individual integrals :math:`\langle\Phi_i|f_{i,j}|\Phi_j\rangle`.
        :type summed: bool, default is ``False``.
        :param component: Request only the i-th component of the result. Remember that :math:`i \in [0, N^2-1]`.
        :param diag_component: Request only the i-th component from the diagonal entries, here :math:`i \in [0, N-1]`.
                               Note that ``component`` takes precedence over ``diag_component`` if both are supplied. (Which is discouraged)
        :return: The value of :math:`\langle\Psi|f|\Psi\rangle`. This is either a scalar value or a list of :math:`N^2` scalar elements.
        """
        eps = packet.get_eps()

        nodes = self.transform_nodes(packet.get_parameters(), eps)
        weights = self._QR.get_weights()

        N = packet.get_number_components()
        K = [ bs.get_basissize() for bs in packet.get_basis_shape() ]

        coeffs = packet.get_coefficients()

        # Avoid unnecessary computations of other components
        if component is not None:
            rows = [ component // N ]
            cols = [ component % N ]
        elif diag_component is not None:
            rows = [ diag_component ]
            cols = [ diag_component ]
        else:
            rows = xrange(N)
            cols = xrange(N)

        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        #       For this, 'operator' must support the component=(r,c) option.
        if operator is None:
            # Operator is None is interpreted as identity transformation
            nn = self._QR.get_number_nodes()
            operator = lambda nodes, component=None: ones((nn,)) if component[0] == component[1] else zeros((nn,))
            values = [ operator(nodes, component=(r,c)) for r in xrange(N) for c in xrange(N) ]
        else:
            values = operator(nodes)

        # Evaluate only the bases we need:
        bases = [ None for n in xrange(N) ]

        for row in rows:
            if bases[row] is None:
                bases[row] = packet.evaluate_basis_at(nodes, component=row)

        for col in cols:
            if bases[col] is None:
                bases[col] = packet.evaluate_basis_at(nodes, component=col)

        # Compute the quadrature
        result = []
        for row in rows:
            for col in cols:
                M = zeros((K[row],K[col]), dtype=complexfloating)

                factor = squeeze(eps * weights * values[row*N + col])

                # Summing up matrices over all quadrature nodes
                for r in xrange(self._QR.get_number_nodes()):
                    M += factor[r] * outer(conjugate(bases[row][:,r]), bases[col][:,r])

                # And include the coefficients as conj(c)*M*c
                result.append(dot(conjugate(coeffs[row]).T, dot(M,coeffs[col])))

        if summed is True:
            result = sum(result)
        elif component is not None:
            # Do not return a list for quadrature of specific single components
            result = result[0]

        return result


    def build_matrix(self, packet, operator=None):
        r"""Calculate the matrix elements of :math:`\langle\Psi|f|\Psi\rangle`
        for a general function :math:`f(x)` with :math:`x \in \mathbb{R}^D`.

        :param packet: The wavepacket :math:`\Psi`.
        :param operator: A real-valued function :math:`f(q, x): \mathbb{R} \times \mathbb{R}^D \rightarrow \mathbb{R}^{N \times N}`.
        :return: A square matrix of size :math:`\sum_i^N |\mathcal{K}_i| \times \sum_j^N |\mathcal{K}_j|`.
        """
        eps = packet.get_eps()

        nodes = self.transform_nodes(packet.get_parameters(), eps)
        weights = self._QR.get_weights()

        N = packet.get_number_components()
        K = [ bs.get_basissize() for bs in packet.get_basis_shape() ]

        bases = [ packet.evaluate_basis_at(nodes, component=n) for n in xrange(N) ]

        # The partition scheme of the block vectors and block matrix
        partition = [0] + list(cumsum(K))

        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        #       For this, 'operator' must support the component=(r,c) option.
        if operator is None:
            # Operator is None is interpreted as identity transformation
            nn = self._QR.get_number_nodes()
            operator = lambda nodes, component=None: ones((nn,)) if component[0] == component[1] else zeros((nn,))
            values = [ operator(nodes, component=(r,c)) for r in xrange(N) for c in xrange(N) ]
        else:
            # TODO: operator should be only f(nodes)
            q, P, Q, P, S = packet.get_parameters()
            values = operator(q, nodes)

        result = zeros((sum(K),sum(K)), dtype=complexfloating)

        for row in xrange(N):
            for col in xrange(N):
                factor = squeeze(eps * weights * values[row*N + col])

                M = zeros((K[row],K[col]), dtype=complexfloating)

                # Sum up matrices over all quadrature nodes and
                # remember to slice the evaluated basis appropriately
                for r in xrange(self._QR.get_number_nodes()):
                    M += factor[r] * outer(conjugate(bases[row][:,r]), bases[col][:,r])

                result[partition[row]:partition[row+1], partition[col]:partition[col+1]] = M

        return result
