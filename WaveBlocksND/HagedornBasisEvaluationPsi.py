"""The WaveBlocks Project

The basic common algorithms for evaluation Hagedorn basis functions
of the new kind.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import complexfloating, dot, vstack, zeros, identity, diag, real
from numpy.linalg import eigh
from scipy import sqrt
from scipy.linalg import det, polar

from WaveBlocksND.HagedornBasisEvaluationCommon import HagedornBasisEvaluationCommon

__all__ = ["HagedornBasisEvaluationPsi"]


class HagedornBasisEvaluationPsi(HagedornBasisEvaluationCommon):
    r"""
    """

    def evaluate_basis_at(self, grid, component, *, prefactor=False):
        r"""Evaluate the basis functions :math:`\psi_k` recursively at the given nodes :math:`\gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:meth:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: Boolean, default is ``False``.
        :return: A two-dimensional ndarray :math:`H` of shape :math:`(|\mathfrak{K}_i|, |\Gamma|)` where
                 the entry :math:`H[\mu(k), i]` is the value of :math:`\psi_k(\gamma_i)`.
        """
        D = self._dimension

        bas = self._basis_shapes[component]
        bs = self._basis_sizes[component]

        # The grid
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes()
        nn = grid.get_number_nodes(overall=True)

        # Allocate the storage array
        psi = zeros((bs, nn), dtype=complexfloating)

        # Precompute some constants
        Pi = self.get_parameters(component=component)
        q, p, Q, P, _ = Pi

        # Transformation to {w} basis
        _, PA = polar(Q, side='left')
        EW, EV = eigh(real(PA))
        Qinv = dot(diag(1.0 / EW), EV.T)
        QQ = identity(D)

        # Compute the ground state psi_0 via direct evaluation
        mu0 = bas[tuple(D * [0])]
        psi[mu0, :] = self._evaluate_psi0(component, nodes, prefactor=False)

        # Compute all higher order states psi_k via recursion
        for d in range(D):
            # Iterator for all valid index vectors k
            indices = bas.get_node_iterator(mode="chain", direction=d)

            for k in indices:
                # Current index vector
                ki = vstack(k)

                # Access predecessors
                psim = zeros((D, nn), dtype=complexfloating)

                for j, kpj in bas.get_neighbours(k, selection="backward"):
                    mukpj = bas[kpj]
                    psim[j, :] = psi[mukpj, :]

                # Compute 3-term recursion
                p1 = (nodes - q) * psi[bas[k], :]
                p2 = sqrt(ki) * psim

                t1 = sqrt(2.0 / self._eps**2) * dot(Qinv[d, :], p1)
                t2 = dot(QQ[d, :], p2)

                # Find multi-index where to store the result
                kped = bas.get_neighbours(k, selection="forward", direction=d)

                # Did we find this k?
                if len(kped) > 0:
                    kped = kped[0]

                    # Store computed value
                    psi[bas[kped[1]], :] = (t1 - t2) / sqrt(ki[d] + 1.0)

        if prefactor is True:
            psi = psi / self._get_sqrt(component)(det(Q))

        return psi


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

        # Transformation to {w} basis
        _, PA = polar(Q, side='left')
        EW, EV = eigh(real(PA))
        Qinv = dot(diag(1.0 / EW), EV.T)
        QQ = identity(D)

        # The basis shape
        bas = self._basis_shapes[component]
        Z = tuple(D * [0])

        # Book keeping
        todo = []
        newtodo = [Z]
        olddelete = []
        delete = []
        tmp = {}

        # The grid nodes
        grid = self._grid_wrap(grid)
        nn = grid.get_number_nodes(overall=True)
        nodes = grid.get_nodes()

        # Evaluate psi0
        tmp[Z] = self._evaluate_psi0(component, nodes, prefactor=False)
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

                # Access predecessors
                psim = zeros((D, nn), dtype=complexfloating)
                for j, kpj in bas.get_neighbours(k, selection="backward"):
                    psim[j, :] = tmp[kpj]

                # Compute the neighbours
                for d, n in bas.get_neighbours(k, selection="forward"):
                    if n not in tmp.keys():
                        # Compute 3-term recursion
                        p1 = (nodes - q) * tmp[k]
                        p2 = sqrt(ki) * psim

                        t1 = sqrt(2.0 / self._eps**2) * dot(Qinv[d, :], p1)
                        t2 = dot(QQ[d, :], p2)

                        # Store computed value
                        tmp[n] = (t1 - t2) / sqrt(ki[d] + 1.0)
                        # And update the result
                        psi = psi + self._coefficients[component][bas[n], 0] * tmp[n]

                        newtodo.append(n)
                delete.append(k)

        if prefactor is True:
            psi = psi / self._get_sqrt(component)(det(Q))

        return psi
