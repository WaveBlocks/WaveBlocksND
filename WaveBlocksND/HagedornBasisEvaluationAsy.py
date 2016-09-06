"""The WaveBlocks Project

The basic common algorithms for evaluation Hagedorn basis functions
of the old kind.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import complexfloating, dot, zeros, conjugate, exp, real, imag, array
from scipy import sqrt
from scipy.linalg import det, inv

from WaveBlocksND.HagedornBasisEvaluationCommon import HagedornBasisEvaluationCommon
from WaveBlocksND.hermite import hermite_asy


__all__ = ["HagedornBasisEvaluationAsy"]


class HagedornBasisEvaluationAsy(HagedornBasisEvaluationCommon):
    r"""
    """

    def evaluate_basis_at(self, grid, component, *, prefactor=False):
        r"""Evaluate the basis functions :math:`\phi_k` recursively at the given nodes :math:`\gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :return: A two-dimensional ndarray :math:`H` of shape :math:`(|\mathfrak{K}_i|, |\Gamma|)` where
                 the entry :math:`H[\mu(k), i]` is the value of :math:`\phi_k(\gamma_i)`.
        """
        D = self._dimension

        assert D == 1

        bas = self._basis_shapes[component]
        bs = self._basis_sizes[component]

        # The grid
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes()
        nn = grid.get_number_nodes(overall=True)

        # Precompute some constants
        eps = self._eps

        Pi = self.get_parameters(component=component)
        q, p, Q, P, _ = Pi

        Qinv = inv(Q)
        Qbar = conjugate(Q)
        Qabs = dot(Q, Qbar)

        # Transform nodes
        Q0 = inv(dot(Q, conjugate(Q.T)))
        nodest = dot(Q0, nodes - q) / eps
        nodest = real(nodest)

        # Compute R(y)
        Ry = (exp(0.5j * nodest**2 * (real(P) * real(Q) + imag(P) * imag(Q))) *
              exp(1.0j / eps * p * Qabs * nodest))

        # Allocate the storage array
        phi = zeros((bs, nn), dtype=complexfloating)

        # Compute all states phi_k(x) via an asymptotic expansion of h_k(x)
        for k in bas.get_node_iterator(mode="lex"):
            # Current index vector
            ki = array(k)[0]

            # Attention: Q, Qbar are scalars here
            Qk = Qinv**(ki / 2.0) * Qbar**(ki / 2)

            # Hermite function evaluation
            hk = hermite_asy(ki, nodest)

            # Assign
            phi[bas[k], :] = Qk * hk * Ry / sqrt(eps)

        if prefactor is True:
            phi = phi / self._get_sqrt(component)(det(Q))

        return phi


    def slim_recursion(self, grid, component, *, prefactor=False):
        r"""Evaluate the Hagedorn wavepacket :math:`\Psi` at the given nodes :math:`\gamma`.
        This routine is a slim version compared to the full basis evaluation. At every moment
        we store only the data we really need to compute the next step until we hit the highest
        order basis functions.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :return: A list of arrays or a single array containing the values of the :math:`\Phi_i`
                 at the nodes :math:`\gamma`.

        Note that this function does not include the global phase :math:`\exp(\frac{i S}{\varepsilon^2})`.
        """
        D = self._dimension

        assert D == 1

        bas = self._basis_shapes[component]

        # The grid
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes()
        nn = grid.get_number_nodes(overall=True)

        # Precompute some constants
        eps = self._eps

        Pi = self.get_parameters(component=component)
        q, p, Q, P, _ = Pi

        Qinv = inv(Q)
        Qbar = conjugate(Q)
        Qabs = dot(Q, Qbar)

        # Transform nodes
        Q0 = inv(dot(Q, conjugate(Q.T)))
        nodest = dot(Q0, nodes - q) / eps
        nodest = real(nodest)

        # Compute R(y)
        Ry = (exp(0.5j * nodest**2 * (real(P) * real(Q) + imag(P) * imag(Q))) *
              exp(1.0j / eps * p * Qabs * nodest))

        # Allocate the storage array
        phi = zeros((nn,), dtype=complexfloating)

        # Compute all states phi_k(x) via an asymptotic expansion of h_k(x)
        for k in bas.get_node_iterator(mode="lex"):
            # Current index vector
            ki = array(k)[0]

            # Attention: Q, Qbar are scalars here
            Qk = Qinv**(ki / 2.0) * Qbar**(ki / 2)

            # Hermite function evaluation
            hk = hermite_asy(ki, nodest)

            # Assign
            ck = self._coefficients[component][bas[k], 0]
            phi = phi + ck * Qk * hk * Ry / sqrt(eps)

        if prefactor is True:
            phi = phi / self._get_sqrt(component)(det(Q))

        return phi
