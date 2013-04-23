"""The WaveBlocks Project

Use a symbolic exact formula for computing the inner product
between two semi-classical wavepackets. The formula is
constructed explicitely for the inhomogeneous case.

@author: R. Bourquin
@copyright: Copyright (C) 2013 R. Bourquin
@license: Modified BSD License
"""

from numpy import squeeze, conjugate, sqrt, ones, zeros, complexfloating
from scipy import exp
from scipy.misc import factorial
from scipy.special import binom
from scipy.special.orthogonal import eval_hermite

from Quadrature import Quadrature

__all__ = ["SymbolicIntegral"]


class SymbolicIntegral(Quadrature):
    r"""
    """

    def __init__(self, *unused, **kunused):
        r"""
        """
        # Drop any argument, we do not need a qr instance.


    def __str__(self):
        return "Inhomogeneous inner product computed using a symbolic exact formula."


    def get_description(self):
        r"""Return a description of this integral object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "SymbolicIntegral"
        return d


    def initialize_packet(self, pacbra, packet=None):
        r"""Provide the wavepacket parts of the inner product to evaluate.
        Since the formula is for the inhomogeneous case explicitly, different
        wavepackets can be used for the 'bra' as well as the 'ket' part.

        :param pacbra: The packet that is used for the 'bra' part.
        :param packet: The packet that is used for the 'ket' part.
        :raises: :py:class:`ValueError` if the dimension of :math:`\Psi` is not 1.
        """
        # Allow to ommit the ket if it is the same as the bra
        if packet is None:
            packet = pacbra

        if not pacbra.get_dimension() == 1:
            raise ValueError("The 'SymbolicIntegral' applies in the 1D case only.")

        self._pacbra = pacbra
        self._packet = packet


    def initialize_operator(self, operator=None, matrix=False):
        r"""Provide the operator part of the inner product to evaluate.
        This function initializes the operator used for quadratures
        and for building matrices.

        Note that the symbolic solution can not handle operators at all.

        :param operator: The operator of the inner product.
                         If `None` a suitable identity is used.
        :param matrix: Set this to ``True`` (Default is ``False``) in case
                       we want to compute the matrix elements.
                       For nasty technical reasons we can not yet unify
                       the operator call syntax.
        """
        # Operator is None is interpreted as identity transformation
        if operator is None:
            self._operator = lambda nodes, dummy, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
        else:
            raise ValueError("The 'SymbolicIntegral' can not handle operators.")


    def prepare(self, rows, cols):
        r"""Precompute some values needed for evaluating the integral
        :math:`\langle \Phi_i | \Phi^\prime_j \rangle` or the corresponding
        matrix over the basis functions of :math:`\Phi_i` and :math:`\Phi^\prime_j`.
        Note that this function does nothing in the current implementation.

        :param rows: A list of all :math:`i` with :math:`0 \leq i \leq N`
                     selecting the :math:`\Phi_i` for which te precompute values.
        :param cols: A list of all :math:`j` with :math:`0 \leq j \leq N`
                     selecting the :math:`\Phi^\prime_j` for which te precompute values.
        """
        pass


    def exact_result_ground(self, Pibra, Piket, eps):
        q1, p1, Q1, P1 = Pibra
        q2, p2, Q2, P2 = Piket
        hbar = eps**2
        X = Q2*conjugate(P1) - P2*conjugate(Q1)
        I = sqrt(-2.0j/X) * exp( 1.0j/(2*hbar) * (Q2*conjugate(Q1)*(p2 - p1)**2 + P2*conjugate(P1)*(q2 - q1)**2) / X
                                -1.0j/hbar *     ((q2 - q1)*(Q2*conjugate(P1)*p2 - P2*conjugate(Q1)*p1)) / X
                               )
        return I


    def exact_result_higher(self, Pibra, Piket, eps, l, k):
        # < phi_l[Pi1] | phi_k[Pi2] >
        q1, p1, Q1, P1 = Pibra
        q2, p2, Q2, P2 = Piket

        I0 = self.exact_result_ground(Pibra, Piket, eps)
        pf = (1.0/sqrt(factorial(k)*factorial(l)) * 2**(-(l+k)/2.0) * I0 *
              (1.0j*conjugate(P1)*Q2-1.0j*conjugate(Q1)*P2)**(-(l+k)/2.0))

        # TODO: Note: formula currently fails for non-inhomogeneous case
        #       because of divisions by zero in the two args below.
        arg1 = ((1.0j*conjugate(P1)*(q1-q2) - 1.0j*conjugate(Q1)*(p1-p2)) /
                (sqrt(1.0j*conjugate(Q2*P1) - 1.0j*conjugate(Q1*P2)) *
                 sqrt(1.0j*conjugate(P1)*Q2 - 1.0j*conjugate(Q1)*P2)))

        arg2 = ((1.0j*Q2*(p1-p2) - 1.0j*P2*(q1 - q2)) /
                (sqrt(1.0j*Q2*P1 - 1.0j*Q1*P2)*sqrt(1.0j*conjugate(P1)*Q2 - 1.0j*conjugate(Q1)*P2)))

        S = 0.0j

        for j in xrange(0, min(l,k)+1):
            s = (binom(l,j)*binom(k,j) * factorial(j) * 4**j *
                 (1.0j*conjugate(Q2*P1) - 1.0j*conjugate(Q1*P2))**((k-j)/2.0) *
                 (1.0j*Q2*P1 - 1.0j*Q1*P2)**((l-j)/2.0))

            h1 = eval_hermite(k-j,  1.0/eps * arg1)
            h2 = eval_hermite(l-j, -1.0/eps * arg2)

            S = S + s*h1*h2
        #else:
        #    S = 1.0

        Ikl = pf * S
        return squeeze(Ikl)


    def perform_quadrature(self, row, col):
        r"""Evaluates the integral :math:`\langle \Phi_i | \Phi^\prime_j \rangle`
        by an exact symbolic formula.

        :param row: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`.
        :param row: The index :math:`j` of the component :math:`\Phi^\prime_j` of :math:`\Psi^\prime`.
        :return: A single complex floating point number.
        """
        eps = self._packet.get_eps()

        Pibra = self._pacbra.get_parameters(component=row)
        Piket = self._packet.get_parameters(component=col)
        cbra = self._pacbra.get_coefficient_vector(component=row)
        cket = self._packet.get_coefficient_vector(component=col)
        Kbra = self._pacbra.get_basis_shapes(component=row)
        Kket = self._packet.get_basis_shapes(component=col)

        result = 0.0j

        for r in Kbra:
            for c in Kket:
                cr = cbra[Kbra[r],0]
                cc = cket[Kket[c],0]
                i = self.exact_result_higher(Pibra[:4], Piket[:4], eps, r[0], c[0])
                result = result + conjugate(cr) * cc * i

        phase = exp(1.0j/eps**2 * (Piket[4]-conjugate(Pibra[4])))
        return phase * result


    def perform_build_matrix(self, row, col):
        r"""Computes the matrix elements :math:`\langle\Phi_i |\Phi^\prime_j\rangle`
        by an exact symbolic formula.

        :param row: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`.
        :param row: The index :math:`j` of the component :math:`\Phi^\prime_j` of :math:`\Psi^\prime`.
        :return: A complex valued matrix of shape :math:`|\mathcal{K}_i| \times |\mathcal{K}^\prime_j|`.
        """
        eps = self._packet.get_eps()

        Pibra = self._pacbra.get_parameters(component=row)
        Piket = self._packet.get_parameters(component=col)
        Kbra = self._pacbra.get_basis_shapes(component=row)
        Kket = self._packet.get_basis_shapes(component=col)

        M = zeros((Kbra.get_basis_size(),Kket.get_basis_size()), dtype=complexfloating)

        for r in Kbra:
            for c in Kket:
                M[Kbra[r], Kket[c]] = self.exact_result_higher(Pibra[:4], Piket[:4], eps, r[0], c[0])

        phase = exp(1.0j/eps**2 * (Piket[4]-conjugate(Pibra[4])))
        return phase * M
