"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2013 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, ones, complexfloating, squeeze, conjugate, dot, outer
from scipy.linalg import sqrtm, inv #, svd, diagsvd

from Quadrature import Quadrature

__all__ = ["DirectHomogeneousQuadrature"]


class DirectHomogeneousQuadrature(Quadrature):
    r"""
    """

    def __init__(self, QR=None):
        # Pure convenience to allow setting of quadrature rule in constructor
        if QR is not None:
            self.set_qr(QR)
        else:
            self._QR = None


    def initialize_packet(self, packet):
        r"""
        """
        self._packet = packet

        # Adapt the quadrature nodes and weights
        eps = self._packet.get_eps()
        self._nodes = self.transform_nodes(self._packet.get_parameters(), eps)
        self._weights = self._QR.get_weights()

        # Force a call of 'preprare'
        self._values = None
        self._bases = None
        self._coeffs = None


    def initialize_operator(self, operator=None):
        r"""
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        #       For this, 'operator' must support the 'component=(r,c)' option.
        N  = self._packet.get_number_components()
        if operator is None:
            # Operator is None is interpreted as identity transformation
            self._operator = lambda nodes, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
            self._values = tuple([ self._operator(self._nodes, entry=(r,c)) for r in xrange(N) for c in xrange(N) ])
        else:
            self._operator = operator
            self._values = tuple( self._operator(self._nodes) )

        # Recheck what we got
        assert type(self._values) is tuple
        assert len(self._values) == N**2


    def initialize_operator_matrix(self, operator=None):
        r"""
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        # For this, 'operator' must support the 'entry=(r,c)' option.
        N  = self._packet.get_number_components()
        if operator is None:
            # Operator is None is interpreted as identity transformation
            self._operator = lambda nodes, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
            self._values = tuple([ self._operator(self._nodes, entry=(r,c)) for r in xrange(N) for c in xrange(N) ])
        else:
            self._operator = operator
            # TODO: operator should be only f(nodes) but we can not fix this currently
            q, p, Q, P, S = self._packet.get_parameters()
            self._values = tuple( self._operator(self._nodes, q) )


    def prepare(self, rows, cols):
        r"""
        """
        # Evaluate only the bases we need
        N  = self._packet.get_number_components()
        bases = [ None for n in xrange(N) ]

        for row in rows:
            if bases[row] is None:
                bases[row] = self._packet.evaluate_basis_at(self._nodes, component=row, prefactor=False)

        for col in cols:
            if bases[col] is None:
                bases[col] = self._packet.evaluate_basis_at(self._nodes, component=col, prefactor=False)

        self._bases = bases

        # Coefficients
        self._coeffs = self._packet.get_coefficients()


    def prepare_for_row(self, row, col):
        pass

    def preprare_for_col(self, row, col):
        pass


    def transform_nodes(self, Pi, eps, QR=None):
        r"""Transform the quadrature nodes :math:`\gamma` such that they
        fit the given wavepacket :math:`\Phi\left[\Pi\right]`.

        :param Pi: The parameter set :math:`\Pi` of the wavepacket.
        :param eps: The value of :math:`\varepsilon` of the wavepacket.
        :param QR: An optional quadrature rule :math:`\Gamma = (\gamma, \omega)` providing
                   the nodes. If not given the internal quadrature rule will be used.
        :return: A two-dimensional ndarray of shape :math:`(D, |\Gamma|)` where
                 :math:`|\Gamma|` denotes the total number of quadrature nodes.
        """
        if QR is None:
            QR = self._QR

        if QR["transform"] is not True:
            return QR.get_nodes()

        q, p, Q, P, S = Pi

        # Compute the affine transformation
        Q0 = inv(dot(Q, conjugate(Q.T)))
        Qs = inv(sqrtm(Q0))

        # TODO: Avoid sqrtm and inverse computation, use svd
        # Untested code:
        # D = Q.shape[0]
        # U, S, Vh = svd(inv(Q))
        # Sinv = 1.0 / S
        # Sinv = diagsvd(Sinv, D, D)
        # Qs = dot(conjugate(Vh.T), dot(Sinv, Vh))

        # And transform the nodes
        nodes = q + eps * dot(Qs, QR.get_nodes())
        return nodes.copy()


    def perform_quadrature(self, row, col):
        r"""
        """
        D = self._packet.get_dimension()
        eps = self._packet.get_eps()

        N  = self._packet.get_number_components()
        Kr = self._packet.get_basis_shapes(component=row).get_basis_size()
        Kc = self._packet.get_basis_shapes(component=col).get_basis_size()

        M = zeros((Kr,Kc), dtype=complexfloating)

        factor = squeeze(eps**D * self._weights * self._values[row*N + col])

        # Summing up matrices over all quadrature nodes
        for r in xrange(self._QR.get_number_nodes()):
            M += factor[r] * outer(conjugate(self._bases[row][:,r]), self._bases[col][:,r])

        # And include the coefficients as conj(c).T*M*c
        return dot(conjugate(self._coeffs[row]).T, dot(M, self._coeffs[col]))


    def perform_build_matrix(self, row, col):
        r"""
        """
        D = self._packet.get_dimension()
        eps = self._packet.get_eps()

        N  = self._packet.get_number_components()
        Kr = self._packet.get_basis_shapes(component=row).get_basis_size()
        Kc = self._packet.get_basis_shapes(component=col).get_basis_size()

        M = zeros((Kr,Kc), dtype=complexfloating)

        factor = squeeze(eps**D * self._weights * self._values[row*N + col])

        # Sum up matrices over all quadrature nodes
        for r in xrange(self._QR.get_number_nodes()):
            M += factor[r] * outer(conjugate(self._bases[row][:,r]), self._bases[col][:,r])

        return M
