"""The WaveBlocks Project

Use a symbolic exact formula for computing the inner product
between two semi-classical wavepackets. The formula is
constructed explicitly for the inhomogeneous case.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import squeeze, conjugate, sqrt, arange, isnan
from scipy.misc import factorial

from InnerProduct import InnerProductException
from SymbolicIntegral import SymbolicIntegral

__all__ = ["SymbolicIntegral0"]


class SymbolicIntegral0(SymbolicIntegral):
    r"""
    """

    def __str__(self):
        return "Inhomogeneous inner product <0|l> computed using a symbolic exact formula."


    def get_description(self):
        r"""Return a description of this integral object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "SymbolicIntegral0"
        return d


    def exact_result_higher(self, Pibra, Piket, eps, k, l):
        r"""Compute the overlap integral :math:`\langle \phi_0 | \phi_l \rangle` of
        two states :math:`\phi_0` and :math:`\phi_l` by using the symbolic formula:

        .. math::
            \langle \phi_0 | \phi_l \rangle =
            \frac{1}{\sqrt{l!}} 2^{-\frac{l}{2}} \langle \phi_0 | \phi_0 \rangle \cdot
            \left(i \overline{ P_1} Q_2 - i \overline{Q_1} P_2\right)^{-\frac{l}{2}} \cdot
            \left(i \overline{Q_2 P_1} - i\overline{Q_1 P_2}\right)^{\frac{l}{2}}
            \cdot
            \\
            H_{l}\left(
              \frac{1}{\varepsilon}
              \frac{\overline{ P_1}\left(q_1-q_2\right)-\overline{Q_1}\left(p_1-p_2\right)}
                   {\sqrt{\overline{Q_2 P_1}-\overline{Q_1 P_2}}\sqrt{\overline{ P_1}Q_2-\overline{Q_1} P_2}}\right)

        Note that this is an internal method and usually there is no
        reason to call it from outside.

        :param Pibra: The parameter set :math:`\Pi = \{q_1,p_1,Q_1,P_1\}` of the bra :math:`\langle \phi_0 |`.
        :param Piket: The parameter set :math:`\Pi^\prime = \{q_2,p_2,Q_2,P_2\}` of the ket :math:`| \phi_l \rangle`.
        :param eps: The semi-classical scaling parameter :math:`\varepsilon`.
        :param k: Index :math:`k=0` of the wavepacket basis function :math:`\phi_0`.
        :param l: Index :math:`l` of the wavepacket basis function :math:`\phi_l`.
        :return: The value of the integral :math:`\langle \phi_0 | \phi_l \rangle`.
        """
        # Do not raise because wavepackets might have a larger basis shape
        # with only the coefficient 0 having a non-zero value. Returning
        # a value 0.0j for k > 0 seems to be a safe default behaviour.
        if not k == 0:
            return 0.0j

        # < phi_0[Pi1] | phi_l[Pi2] >
        q1, p1, Q1, P1 = Pibra
        q2, p2, Q2, P2 = Piket

        # If both parameter sets are identical, we are back in the
        # homogeneous case where the phi are orthonormal.
        # This early returns just serves to avoid introducing NaN
        # elements in arg1 and arg2 further below. It is allowed to
        # compare floats on exact equality because the inhomogeneous
        # formula works even if the floats differ only by tiny amounts.
        if q1 == q2 and p1 == p2 and Q1 == Q2 and P1 == P2:
            return 1.0 if l == 0 else 0.0

        # TODO: Note that the formula can still fail if Q1 = Q2 and P1 = P2
        #       but q1 \neq q2 and p1 \neq p2.
        pf = (self._f[0,l] * 2**(-l/2.0) * self._I0 *
              (1.0j*conjugate(P1)*Q2-1.0j*conjugate(Q1)*P2)**(-l/2.0))

        S = self._pfl[l] * self._Hl[l]

        Ikl = pf * S
        return squeeze(Ikl)


    def _cache_factors(self, Pibra, Piket, Kbra, Kket, eps):
        r"""Cache some summands to speed up the computation of the sum.

        :param Pibra: The parameter set :math:`\Pi` of the bra :math:`\langle \Phi |`.
        :param Piket: The parameter set :math:`\Pi^\prime` of the ket :math:`| \Phi^\prime \rangle`.
        :param Kbra: The basis shape :math:`\mathfrak{K}` of the bra :math:`\langle \Phi |`.
        :type Kbra: A :py:class:`BasisShape` instance.
        :param Kket: The basis shape :math:`\mathfrak{K}^\prime` of the ket :math:`| \Phi^\prime \rangle`.
        :type Kket: A :py:class:`BasisShape` instance.
        :param eps: The semi-classical scaling parameter :math:`\varepsilon`.
        """
        q1, p1, Q1, P1 = Pibra
        q2, p2, Q2, P2 = Piket

        # If both parameter sets are identical, we are back in the homogeneous case.
        if q1 == q2 and p1 == p2 and Q1 == Q2 and P1 == P2:
            self._Hl = None

        # We have k in [0, 1, ..., K-1] where |K| is the basis size
        # hence K-1 is the maximal index.
        L = Kket.get_basis_size()

        makl = L

        # Factorials
        f = factorial(arange(makl))
        self._f = 1.0 / sqrt(f[:L].reshape(1,-1))

        # Note: formula currently fails for non-inhomogeneous case
        #       because of divisions by zero in the two args below.
        argl = ((1.0j*conjugate(P1)*(q1-q2) - 1.0j*conjugate(Q1)*(p1-p2)) /
                (sqrt(1.0j*conjugate(Q2*P1) - 1.0j*conjugate(Q1*P2)) *
                 sqrt(1.0j*conjugate(P1)*Q2 - 1.0j*conjugate(Q1)*P2)))

        # TODO: Better test for failure?
        if self._doraise and isnan(squeeze(argl)):
            raise InnerProductException("Symbolic formula failed due to Q_k = Q_l and P_k = P_l.")

        # The parameter j varies in the range [0, 1, ..., min(K-1,L-1)]
        # hence we have that k-j can be in [K-1, K-2, ..., K-1-min(K-1,L-1)]
        # and similar for l-j we have [L-1, L-2, ..., L-1-min(K-1,L-1)]
        # where both K-1-min(K-1,L-1) and L-1-min(K-1,L-1) are non-negative.
        self._Hl = self._evaluate_hermite(L-1, 1.0/eps * argl)

        il = arange(L).reshape(-1,1)
        self._pfl = ((1.0j*conjugate(Q2*P1) - 1.0j*conjugate(Q1*P2)) ** (il/2.0)).reshape(L)

        # And the groundstate value
        self._I0 = self.exact_result_ground(Pibra, Piket, eps)
