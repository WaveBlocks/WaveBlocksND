"""The WaveBlocks Project

The basic common algorithms for evaluation Hagedorn basis functions
of the old kind.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import complexfloating, dot, vstack, zeros, conjugate, exp, real, imag, squeeze, array
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

        # Allocate the storage array
        phi = zeros((bs, nn), dtype=complexfloating)

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

        # Compute all higher order states phi_k via asymptotic expansion
        for d in range(D):
            # Iterator for all valid index vectors k
            indices = list(bas.get_node_iterator(mode="chain", direction=d))
            print(list(indices))

            for k in indices:
                # Current index vector
                ki = array(k)
                print(ki)

                # Attention: Q, Qbar are scalars here
                Qk = dot(Qinv**((ki + 1) / 2.0), Qbar**(ki / 2))

                # Compute R(y)
                Rk = (exp(0.5j * nodest**2 * (real(P) * real(Q) + imag(P) * imag(Q))) *
                      exp(1.0j / eps * p * Qabs * nodest))

                # Hermite function evaluation
                hk = hermite_asy(squeeze(ki), nodest)

                # Assign
                value = 1.0 / sqrt(eps) * Qk * hk * Rk
                phi[bas[k], :] = value

        # if prefactor is True:
        #    phi = phi / self._get_sqrt(component)(det(Q))

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

        # Precompute some constants
        Pi = self.get_parameters(component=component)
        q, p, Q, P, _ = Pi

        Qinv = inv(Q)
        Qbar = conjugate(Q)
        QQ = dot(Qinv, Qbar)

        # The basis shape
        bas = self._basis_shapes[component]
        Z = tuple(D * [0])


        # The grid nodes
        grid = self._grid_wrap(grid)
        nn = grid.get_number_nodes(overall=True)
        nodes = grid.get_nodes()

        # Evaluate phi0
        tmp[Z] = self._evaluate_phi0(component, nodes, prefactor=False)
        psi = self._coefficients[component][bas[Z], 0] * tmp[Z]

        # Iterate for higher order states
        while len(newtodo) != 0:
            # Delete results that never will be used again
            for d in olddelete:
                del tmp[d]

            # Exchange queues
            todo = newtodo
            newtodo = []
            olddelete = delete
            delete = []

            # Compute new results
            for k in todo:
                # Center stencil at node k
                ki = vstack(k)


        if prefactor is True:
            psi = psi / self._get_sqrt(component)(det(Q))

        return psi
