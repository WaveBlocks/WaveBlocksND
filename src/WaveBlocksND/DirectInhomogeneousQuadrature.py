"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2013 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, ones, complexfloating, squeeze, imag, conjugate, outer, dot, ndarray
from scipy import exp
from scipy.linalg import sqrtm, inv, det

from Quadrature import Quadrature

__all__ = ["DirectInhomogeneousQuadrature"]


class DirectInhomogeneousQuadrature(Quadrature):
    r"""
    """

    def __init__(self, QR=None):
        # Pure convenience to allow setting of quadrature rule in constructor
        if QR is not None:
            self.set_qr(QR)
        else:
            self._QR = None


    def initialize_packet(self, pacbra, packet=None):
        r"""
        """
        # Allow to ommit the ket if it is the same as the bra
        if packet is None:
            packet = pacbra

        self._pacbra = pacbra
        self._packet = packet

        # Adapt the quadrature nodes and weights
        self._weights = self._QR.get_weights()

        # Force a call of 'preprare'
        self._coeffbra = None
        self._coeffket = None


    def initialize_operator(self, operator=None):
        r"""
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        #       For this, 'operator' must support the 'component=(r,c)' option.
        # Operator is None is interpreted as identity transformation
        if operator is None:
            self._operator = lambda nodes, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
        else:
            self._operator = operator


    def initialize_operator_matrix(self, operator=None):
        r"""
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        # For this, 'operator' must support the 'entry=(r,c)' option.
        # Operator is None is interpreted as identity transformation
        if operator is None:
            self._operatorm = lambda nodes, dummy, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
        else:
            self._operatorm = operator


    def prepare(self, rows, cols):
        r"""
        """
        # Coefficients
        self._coeffbra = self._pacbra.get_coefficients()
        self._coeffket = self._packet.get_coefficients()


    def prepare_for_row(self, row, col):
        pass

    def preprare_for_col(self, row, col):
        pass


    def mix_parameters(self, Pibra, Piket):
        r"""Mix the two parameter sets :math:`\Pi_i` and :math:`\Pi_j`
        from the bra and the ket wavepackets :math:`\Phi\left[\Pi_i\right]`
        and :math:`\Phi^\prime\left[\Pi_j\right]`.

        :param Pibra: The parameter set :math:`\Pi_i` from the bra part wavepacket.
        :param Piket: The parameter set :math:`\Pi_j` from the ket part wavepacket.
        :return: The mixed parameters :math:`q_0` and :math:`Q_S`. (See the theory for details.)
        """
        # <Pibra | ... | Piket>
        (qr, pr, Qr, Pr, Sr) = Pibra
        (qc, pc, Qc, Pc, Sc) = Piket

        # Mix the parameters
        Gr = dot(Pr, inv(Qr))
        Gc = dot(Pc, inv(Qc))

        r = imag(Gc - conjugate(Gr.T))
        s = imag(dot(Gc, qc) - dot(conjugate(Gr.T), qr))

        q0 = dot(inv(r), s)
        Q0 = 0.5 * r

        # Here we can not avoid the matrix root by using svd
        Qs = inv(sqrtm(Q0))

        return (q0, Qs)


    def transform_nodes(self, Pibra, Piket, eps, QR=None):
        r"""Transform the quadrature nodes :math:`\gamma` such that they
        fit the given wavepackets :math:`\Phi\left[\Pi_i\right]` and
        :math:`\Phi^\prime\left[\Pi_j\right]` best.

        :param Pibra: The parameter set :math:`\Pi_i` from the bra part wavepacket.
        :param Piket: The parameter set :math:`\Pi_j` from the ket part wavepacket.
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

        # Mix the parameters to compute the affine transformation
        q0, Qs = self.mix_parameters(Pibra, Piket)

        # And transform the nodes
        nodes = q0 + eps * dot(Qs, QR.get_nodes())
        return nodes.copy()


    def perform_quadrature(self, row, col):
        r"""
        """
        D = self._packet.get_dimension()
        eps = self._packet.get_eps()

        # Packets can also have different basis sizes
        Kbra = self._pacbra.get_basis_shapes(component=row).get_basis_size()
        Kket = self._packet.get_basis_shapes(component=col).get_basis_size()

        Pibra = self._pacbra.get_parameters(component=row)
        Piket = self._packet.get_parameters(component=col)

        M = zeros((Kbra,Kket), dtype=complexfloating)

        # Transform nodes and evaluate bases
        nodes = self.transform_nodes(Pibra, Piket, eps)
        basisr = self._pacbra.evaluate_basis_at(nodes, component=row, prefactor=True)
        basisc = self._packet.evaluate_basis_at(nodes, component=col, prefactor=True)

        # Operator should support the component notation for efficiency
        values = self._operator(nodes, entry=(row,col))

        # Recheck what we got
        assert type(values) is ndarray
        assert values.shape == (1,self._QR.get_number_nodes())

        Pimix = self.mix_parameters(Pibra, Piket)
        factor = squeeze(eps**D * values * self._weights * det(Pimix[1]))

        # Summing up matrices over all quadrature nodes
        for r in xrange(self._QR.get_number_nodes()):
            M += factor[r] * outer(conjugate(basisr[:,r]), basisc[:,r])

        # Compute global phase difference
        phase = exp(1.0j/eps**2 * (Piket[4]-conjugate(Pibra[4])))

        # And include the coefficients as conj(c).T*M*c
        return phase * dot(conjugate(self._coeffbra[row]).T, dot(M, self._coeffket[col]))


    def perform_build_matrix(self, row, col):
        r"""
        """
        D = self._packet.get_dimension()
        eps = self._packet.get_eps()

        # Packets can also have different basis sizes
        Kbra = self._pacbra.get_basis_shapes(component=row).get_basis_size()
        Kket = self._packet.get_basis_shapes(component=col).get_basis_size()

        Pibra = self._pacbra.get_parameters(component=row)
        Piket = self._packet.get_parameters(component=col)

        M = zeros((Kbra,Kket), dtype=complexfloating)

        # Transform nodes and evaluate bases
        nodes = self.transform_nodes(Pibra, Piket, eps)
        basisr = self._pacbra.evaluate_basis_at(nodes, component=row, prefactor=True)
        basisc = self._packet.evaluate_basis_at(nodes, component=col, prefactor=True)

        Pimix = self.mix_parameters(Pibra, Piket)
        # Operator should support the component notation for efficiency
        # TODO: operator should be only f(nodes) but we can not fix this currently
        values = self._operatorm(nodes, Pimix[0], entry=(row,col))

        # Recheck what we got
        assert type(values) is ndarray
        assert values.shape == (1,self._QR.get_number_nodes())

        factor = squeeze(eps**D * values * self._weights * det(Pimix[1]))

        # Sum up matrices over all quadrature nodes
        for r in xrange(self._QR.get_number_nodes()):
            M += factor[r] * outer(conjugate(basisr[:,r]), basisc[:,r])

        # Compute global phase difference
        phase = exp(1.0j/eps**2 * (Piket[4]-conjugate(Pibra[4])))

        return phase * M
