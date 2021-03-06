"""The WaveBlocks Project

Implementation of the unitary transformation between
old-kind and new-kind Hagedorn wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from collections import defaultdict

from numpy import (complexfloating, array, hstack, vstack, split, integer,
                   conjugate, zeros, identity, zeros_like, cumsum, sum)
from scipy import sqrt
from scipy.special import binom


from WaveBlocksND.SimplexShape import SimplexShape
from WaveBlocksND.HagedornWavepacket import HagedornWavepacket
from WaveBlocksND.HagedornWavepacketPsi import HagedornWavepacketPsi

from WaveBlocksND.GaussHermiteQR import GaussHermiteQR
from WaveBlocksND.TensorProductQR import TensorProductQR
from WaveBlocksND.DirectHomogeneousQuadrature import DirectHomogeneousQuadrature

from WaveBlocksND.Combinatorics import multinomial, lattice_points_norm

__all__ = ["HagedornWavepacketTransformPhiPsi"]


class HagedornWavepacketTransformPhiPsi(object):
    r"""Implementation of the unitary transformation between old-kind
    Hagedorn wavepackets :math:`\Phi[\Pi]` and new-kind Hagedorn wavepackets
    :math:`\Psi[\Pi]`.
    """

    def __init__(self):
        pass


    def _build_nu(self, D, order):
        r"""Build the list containing the :math:`\nu_i`.

        :param D: The dimension :math:`D` of the wavepacket.
        :param order: The maximal :math:`l_1` norm of any basis index :math:`\underline{k}`.
        """
        ED = identity(D, dtype=integer)
        nu = [array([zeros((D,))], dtype=integer)]
        for i in range(1, order + 1):
            nu.append(vstack([nu[i - 1] + ed for ed in ED]))

        return nu


    def _build_mu(self, D, order):
        r"""Build the list containing the lists :math:`\mu_i` where each
        :math:`\mu_i` contains the :math:`D`-tuples in reverse lexicographical order.

        :param D: The dimension :math:`D` of the wavepacket.
        :param order: The maximal :math:`l_1` norm of any basis index :math:`\underline{k}`.
        """
        return [list(lattice_points_norm(D, i)) for i in range(order + 1)]


    def _adapt_index_order(self, BS, MU):
        r"""Map the arrangement of the multi-indices :math:`\underline{k}` between the
        native order :math:`\underline{k} \in \mathfrak{K}` and the
        reverse lexicographical order sorted by increasing :math:`l_1` norm.

        :param BS: The basis shape :math:`\mathfrak{K}`.
        :param MU: The list of all :math:`\mu_i`.
        """
        mu = [item for mui in MU for item in mui]
        return array([BS[mu[k]] for k in range(BS.get_basis_size())])


    def overlap(self, D, Pi, eps):
        r"""Compute the overlap matrix:

        .. math::
            \mathbf{K}_{r,c} := \langle \psi_{\underline{e_r}}[\Pi] | \phi_{\underline{e_c}}[\Pi] \rangle

        :param D: The dimension :math:`D`.
        :param Pi: The parameter set :math:`\Pi`.
        :param eps: The semiclassical scaling parameter :math:`\varepsilon`.
        """
        ED = list(lattice_points_norm(D, 1))

        BS = SimplexShape(D, 1)

        WPo = HagedornWavepacket(D, 1, eps)
        WPo.set_parameters(Pi)
        WPo.set_basis_shapes([BS])

        WPn = HagedornWavepacketPsi(D, 1, eps)
        WPn.set_parameters(Pi)
        WPn.set_basis_shapes([BS])

        TPG = TensorProductQR(D * [GaussHermiteQR(2)])
        DHQ = DirectHomogeneousQuadrature(TPG)
        G = DHQ.transform_nodes(Pi, eps)
        W = TPG.get_weights()

        phi = WPo.evaluate_basis_at(G, component=0)
        psi = WPn.evaluate_basis_at(G, component=0)

        K = zeros((D, D), dtype=complexfloating)

        for r, er in enumerate(ED):
            for c, ec in enumerate(ED):
                K[r, c] = eps**D * sum(conjugate(psi[BS[er], :]) * phi[BS[ec], :] * W)

        return K


    def _built_lut(self, nu):
        r"""Build the lookup table for simultaneuos multi-acces to arrays by label lists.

        :param nu: The list of all :math:`\nu_i`.
        """
        lut = defaultdict(list)

        for i, k in enumerate(nu):
            lut[tuple(k)].append(i)

        return lut


    def multiply_Pi_v(self, nu, mu, lut, v):
        r"""Multiply the matrix :math:`\mathbf{P_i}` by a vector :math:`\underline{v}`
        from the right. Do not construct the matrix explicitly.

        :param nu: The list :math:`\nu_i`.
        :param mu: The list :math:`\mu_i`.
        :param lut: The lookup table.
        :param v: The vector :math:`\underline{v} \in \mathcal{X}_i`.
        """
        res = zeros(len(mu), dtype=v.dtype)

        for i in range(len(mu)):
            res[i] = sqrt(multinomial(mu[i])) * v[lut[mu[i]][0]]

        return res


    def multiply_PiT_v(self, nu, mu, lut, v):
        r"""Multiply the matrix :math:`\mathbf{P_i}^{\mathrm{T}}` by a vector
        :math:`\underline{v}` from the right. Do not construct the matrix explicitly.

        :param nu: The list :math:`\nu_i`.
        :param mu: The list :math:`\mu_i`.
        :param lut: The lookup table.
        :param v: The vector :math:`\underline{v} \in \mathbb{C}^{n(D,i)}`.
        """
        res = zeros(nu.shape[0], dtype=v.dtype)

        for i in range(len(mu)):
            res[lut[mu[i]]] = 1.0 / sqrt(multinomial(mu[i])) * v[i]

        return res


    def multiply_kronecker_power_v(self, D, A, v, i):
        r"""Multiply the :math:`i`-th Kronecker power of the matrix :math:`\mathbf{A}`
        by a vector :math:`\underline{v}` from the right.

        :param D: The dimension :math:`D`.
        :param A: The matrix :math:`\mathbf{A}`.
        :param v: The vector :math:`\underline{v}`.
        :param i: The non-negative integer Kronecker power exponent.

        .. note:: The matrix has to be square and of size :math:`D^i`.
        """
        assert i >= 0

        # Select from 'data' exactly 'size' items every 'step' with offset 'base'
        select = lambda data, base, size, step: data.reshape(-1, step)[:, base*size:(base+1)*size]

        if i == 0:
            return v

        w = v
        for j in range(i):
            res = zeros_like(w)
            for r in range(D):
                rhs = select(res, r, D**j, D**(j + 1))
                for c in range(D):
                    lhs = select(w, c, D**j, D**(j + 1))
                    rhs += A[r, c] * lhs

            w = res

        return res


    def _multiply_LOP_v(self, coeffs, LOP, NU, MU, D, J):
        r"""Apply the transformation operator matrix to the coefficient vector.

        :param coeffs: The coefficients vector.
        :param LOP: The transformation operator matrix.
        :param NU: The list of all :math:`\nu_i`.
        :param MU: The list of all :math:`\mu_i`.
        :param D: The dimension :math:`D`.
        :param J: The maximal :math:`l_1` norm of any :math:`\underline{k} \in \mathfrak{K}`.
        """
        nDm = lambda D, m: int(binom(D + m - 1, m))

        s = [nDm(D, i) for i in range(J)]
        Ci = split(coeffs, cumsum(s))
        Cit = []

        for i, ci in enumerate(Ci):
            lut = self._built_lut(NU[i])
            cit = self.multiply_PiT_v(NU[i], MU[i], lut, ci)
            cit = self.multiply_kronecker_power_v(D, LOP, cit, i)
            cit = self.multiply_Pi_v(NU[i], MU[i], lut, cit)
            Cit.append(cit)

        return hstack(Cit).reshape(*coeffs.shape)


    def multiply_T_v(self, coeffs, K, NU, MU, D, J):
        r"""Apply the transformation matrix :math:`\mathbf{T}` to the coefficients :math:`\underline{c}`.

        :param coeffs: The coefficients vector :math:`\underline{c}`.
        :param K: The overlap matrix :math:`\mathbf{K}`.
        :param NU: The list of all :math:`\nu_i`.
        :param MU: The list of all :math:`\mu_i`.
        :param D: The dimension :math:`D`.
        :param J: The maximal :math:`l_1` norm of any :math:`\underline{k} \in \mathfrak{K}`.
        """
        return self._multiply_LOP_v(coeffs, K, NU, MU, D, J)


    def multiply_Tinv_v(self, coeffs, K, NU, MU, D, J):
        r"""Apply the transformation matrix :math:`\mathbf{T}^{-1}`
        to the coefficients :math:`\underline{d}`.

        :param coeffs: The coefficients vector :math:`\underline{d}`.
        :param K: The overlap matrix :math:`\mathbf{K}`.
        :param NU: The list of all :math:`\nu_i`.
        :param MU: The list of all :math:`\mu_i`.
        :param D: The dimension :math:`D`.
        :param J: The maximal :math:`l_1` norm of any :math:`\underline{k} \in \mathfrak{K}`.
        """
        Kinv = K.transpose().conjugate()
        return self._multiply_LOP_v(coeffs, Kinv, NU, MU, D, J)


    def _transform(self, HAWPfrom, LOP, HAWPtoConstructor):
        r"""The transformation logic. This function can transform
        in both directions depending on the linear operator argument.

        :param HAWPfrom: The wavepacket to transform.
        :param LOP: The linear transformation operator, either :math:`\mathbf{T}` or :math:`\mathbf{T}^{-1}`.
        :param HAWPtoConstructor: Constructor of the resulting wavepacket type.
        """
        if not HAWPfrom.get_number_components() == 1:
            raise NotImplementedError("Only scalar wavepackets are supported.")

        D = HAWPfrom.get_dimension()
        BS = HAWPfrom.get_basis_shapes(component=0)

        if not isinstance(BS, SimplexShape):
            raise ValueError("The wavepacket does not have a simplex basis shape.")

        J = BS.get_description()['K']

        NU = self._build_nu(D, J)
        MU = self._build_mu(D, J)

        # Compute small overlap matrix
        Pi = HAWPfrom.get_parameters()
        eps = HAWPfrom.get_eps()
        K = self.overlap(D, Pi, eps)

        # Permutation, adapter to basis index order
        I = self._adapt_index_order(BS, MU)
        cfrom = HAWPfrom.get_coefficients(component=0)
        cto = zeros_like(cfrom)
        cto[I] = LOP(cfrom[I], K, NU, MU, D, J)

        # Set up a to wavepacket
        HAWPto = HAWPtoConstructor(D, 1, eps)
        HAWPto.set_parameters(Pi)
        HAWPto.set_basis_shapes(BS, component=0)
        HAWPto.set_coefficients(cto, component=0)

        return HAWPto


    def transform_phi_to_psi(self, HAWPphi):
        r"""Transform a old-kind wavepacket :math:`\Phi[\Pi]` into a new-kind wavepacket :math:`\Psi[\Pi]`.

        :param HAWPphi: The wavepacket :math:`\Phi[\Pi]` to transform.
        :return: A new wavepacket object representing :math:`\Psi[\Pi]`.
        """
        return self._transform(HAWPphi, self.multiply_T_v, HagedornWavepacketPsi)


    def transform_psi_to_phi(self, HAWPpsi):
        r"""Transform a new-kind wavepacket :math:`\Psi[\Pi]` into a old-kind wavepacket :math:`\Phi[\Pi]`.

        :param HAWPpsi: The wavepacket :math:`\Psi[\Pi]` to transform.
        :return: A new wavepacket object representing :math:`\Phi[\Pi]`.
        """
        return self._transform(HAWPpsi, self.multiply_Tinv_v, HagedornWavepacket)
