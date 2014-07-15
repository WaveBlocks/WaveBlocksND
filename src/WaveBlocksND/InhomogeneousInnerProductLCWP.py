"""The WaveBlocks Project

This file contains code for the delegation of the evaluation of inner products
of two linear combinations of wavepackets. The class defined here can compute
brakets, inner products and expectation values and the matrix elements of an
arbitrary operator.

@author: R. Bourquin
@copyright: Copyright (C) 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, conjugate, transpose, dot

from InnerProduct import InnerProduct

__all__ = ["InhomogeneousInnerProductLCWP"]


class InhomogeneousInnerProductLCWP(InnerProduct):

    def __init__(self, delegate=None, oracle=None):
        r"""
        This class computes the inhomogeneous inner product
        :math:`\langle\Upsilon|f|\Upsilon^\prime\rangle` of two linear combinations
        :math:`\Upsilon` and :math:`\Upsilon^\prime`. The delegate inner
        product class used for computing :math:`\langle\Psi|f|\Psi^\prime\rangle`
        has to be of *inhomogeneous* type.

        :param delegate: The delegate inner product.
        :type delegate: A :py:class:`InnerProduct` subclass instance.
        :param oracle: The sparsity oracle to use. If the variable is ``None``
                       no oracle is used and all integrals are computed.
        """
        # Pure convenience to allow setting of quadrature instance in constructor
        self.set_delegate(delegate)
        self.set_oracle(oracle)


    def __str__(self):
        return "Inhomogeneous inner product of linear combinations computed by " + str(self._delegate)


    def get_description(self):
        r"""Return a description of this inner product object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "InhomogeneousInnerProductLCWP"
        d["delegate"] = self._delegate.get_description()
        return d


    def get_oracle(self):
        r"""Return the sparsity oracle in use or ``None``.
        """
        return self._oracle


    def set_oracle(self, new_oracle):
        r"""Set the sparsity oracle.

        :param new_oracle: The new oracle to use. If the variable is ``None``
                           no oracle is used and all integrals are computed.
        """
        if new_oracle is not None:
            self._oracle = new_oracle
            self._obey_oracle = True
        else:
            self._oracle = None
            self._obey_oracle = False


    def quadrature(self, lcbra, lcket=None, operator=None, component=0):
        r"""Delegates the evaluation of :math:`\langle\Upsilon|f|\Upsilon^\prime\rangle` for a general
        function :math:`f(x)` with :math:`x \in \mathbb{R}^D`.

        :param lcbra: The linear combination :math:`\Upsilon` from the bra with :math:`J` summands :math:`\Psi_j`.
        :param lcket: The linear combination :math:`\Upsilon^\prime` from the ket with :math:`J^\prime` summands :math:`\Psi_j^\prime`.
        :param operator: A matrix-valued function :math:`f(x): \mathbb{R}^D \rightarrow \mathbb{R}^{N \times N^\prime}`.
        :param component: The index :math:`i` of the component :math:`\Phi_j` of :math:`\Psi_j`. If set only those
                          components will be taken into account for the computation.
        :type component: Integer or ``None``, default is ``0``.
        :return: The value of :math:`\langle\Upsilon|f|\Upsilon^\prime\rangle`.
        :type: An :py:class:`ndarray`.
        """

        M = self.build_matrix(lcbra, lcket, operator=operator, component=component)
        cbra = lcbra.get_coefficients()
        cket = lcket.get_coefficients()
        return dot(conjugate(transpose(cbra)), dot(M, cket))


    def build_matrix(self, lcbra, lcket=None, operator=None, component=0):
        r"""Delegates the computation of the matrix elements of :math:`\langle\Upsilon|f|\Upsilon^\prime\rangle`
        for a general function :math:`f(x)` with :math:`x \in \mathbb{R}^D`.
        The matrix is computed without including the coefficients :math:`c_j` and :math:`c_j^\prime`.

        :param lcbra: The linear combination :math:`\Upsilon` from the bra with :math:`J` summands :math:`\Psi_j`.
        :param lcket: The linear combination :math:`\Upsilon^\prime` from the ket with :math:`J^\prime` summands :math:`\Psi_j^\prime`.
        :param operator: A matrix-valued function :math:`f(q, x): \mathbb{R} \times \mathbb{R}^D \rightarrow \mathbb{R}^{N \times N^\prime}`.
        :param component: The index :math:`i` of the component :math:`\Phi_j` of :math:`\Psi_j`. If set only those
                          components will be taken into account for the computation.
        :type component: Integer or ``None``, default is ``0``.
        :return: A matrix of size :math:`J \times J^\prime`.
        :type: An :py:class:`ndarray`.
        """
        # Allow to omit the ket if it is the same as the bra
        if lcket is None:
            lcket = lcbra

        Jbra = lcbra.get_number_packets()
        Jket = lcket.get_number_packets()

        M = zeros((Jbra, Jket), dtype=complexfloating)

        for row, pacbra in enumerate(lcbra.get_wavepackets()):
            for col, packet in enumerate(lcket.get_wavepackets()):
                if self._obey_oracle:
                    if self._oracle.is_not_zero(pacbra, packet, component=component):
                        M[row, col] = self._delegate.quadrature(pacbra, packet, operator=operator, diag_component=component)
                else:
                    M[row, col] = self._delegate.quadrature(pacbra, packet, operator=operator, diag_component=component)

        return M
