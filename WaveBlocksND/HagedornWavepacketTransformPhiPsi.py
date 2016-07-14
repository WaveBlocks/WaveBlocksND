"""The WaveBlocks Project


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
from WaveBlocksND.HagedornWavepacketNew import HagedornWavepacketNew

from WaveBlocksND.GaussHermiteQR import GaussHermiteQR
from WaveBlocksND.TensorProductQR import TensorProductQR
from WaveBlocksND.DirectHomogeneousQuadrature import DirectHomogeneousQuadrature

from WaveBlocksND.Combinatorics import multinomial, lattice_points_norm

__all__ = ["HagedornWavepacketTransformPhiPsi"]


class HagedornWavepacketTransformPhiPsi(object):
    r"""
    """
    def __init__(self):
        pass


    def build_nu(self, D, order):
        r"""Build the list containing the :math:`\nu_i`.

        :param D: The dimension :math:`D` of the wavepacket.
        :param order: The maximal :math:`l_1` norm of any basis index :math:`k`.
        """
        # ED = vsplit(identity(D, dtype=integer), D)
        ED = identity(D, dtype=integer)
        nu = [array([zeros((D,))], dtype=integer)]
        for i in range(1, order + 1):
            nu.append(vstack([nu[i - 1] + ed for ed in ED]))

        return nu


    def build_mu(self, D, order):
        r"""Build the list containing the :math:`\mu_i`.

        :param D: The dimension :math:`D` of the wavepacket.
        :param order: The maximal :math:`l_1` norm of any basis index :math:`k`.
        """
        # TODO: Replace by K-nativ ordering
        return [vstack(lattice_points_norm(D, i)) for i in range(order + 1)]


    def mapit(self, BS, MU, c):
        # TODO: Will resolve in identity
        cm = zeros_like(c)
        K = BS.get_basis_size()
        mu = vstack(MU)
        for k in range(K):
            cm[BS[tuple(mu[k, :])]] = c[k]
        return cm


    def mapitinv(self, BS, MU, c):
        # TODO: Will resolve in identity
        cm = zeros_like(c)
        K = BS.get_basis_size()
        mu = vstack(MU)
        for k in range(K):
            cm[k] = c[BS[tuple(mu[k, :])]]
        return cm


    def overlap(self, D, Pi, eps):
        r"""Compute the overlap matrix:

        .. math::
            \mathbf{K}_{r,c} \assign \langle \psi_{\underline{e_r}}[\Pi] | \phi_{\underline{e_c}}[\Pi] \rangle

        :param D: The dimension :math:`D`.
        :param Pi: The parameter set :math:`\Pi`.
        :param eps: The semiclassical scaling parameter :math:`\varepsilon`.
        """
        ED = list(lattice_points_norm(D, 1))

        # TODO: Recheck and be careful with lima order
        K = SimplexShape(D, 1)

        WPo = HagedornWavepacket(D, 1, eps)
        WPo.set_parameters(Pi)
        WPo.set_basis_shapes([K])

        WPn = HagedornWavepacketNew(D, 1, eps)
        WPn.set_parameters(Pi)
        WPn.set_basis_shapes([K])

        TPG = TensorProductQR(D * [GaussHermiteQR(2)])
        DHQ = DirectHomogeneousQuadrature(TPG)
        G = DHQ.transform_nodes(Pi, eps)
        W = TPG.get_weights()

        phi = WPo.evaluate_basis_at(G, component=0)
        psi = WPn.evaluate_basis_at(G, component=0)

        M = zeros((D, D), dtype=complexfloating)

        for r, er in enumerate(ED):
            for c, ec in enumerate(ED):
                M[r, c] = eps**D * sum(conjugate(psi[K[er], :]) * phi[K[ec], :] * W)

        return M


    def build_lut(self, nu):
        r"""Build the lookup table for simultaneuos multi-acces to arrays by label lists.

        :param nu: The list of all :math:`\nu_i`.
        """
        lut = defaultdict(list)

        for i, k in enumerate(nu):
            lut[tuple(k)].append(i)

        # TODO: Maybe convert lists to arrays?
        return lut


    def multiply_Pi_v(self, nu, mu, lut, v):
        r"""Multiply the matrix :math:`P_i` by a underlinetor :math:`v` from the right.

        :param lut: The lookup table.
        :param v: The vector :math:`v \in \mathcal{X}_i`.
        """
        nmu = mu.shape[0]
        res = zeros(nmu, dtype=v.dtype)

        for i in range(nmu):
            mui = tuple(mu[i, :])
            j = lut[mui][0]
            res[i] = sqrt(multinomial(mui)) * v[j]

        return res


    def multiply_PiT_v(self, nu, mu, lut, v):
        r"""Multiply the matrix :math:`P_i^{T}` by a vector :math:`v` from the right.

        :param nu: The list :math:`\nu`.
        :param mu: The list :math:`\mu`.
        :param lut: The lookup table.
        :param v: The vector :math:`v \in \mathbb{C}^{n(D,i)}`.
        """
        nmu = mu.shape[0]
        nnu = nu.shape[0]
        res = zeros(nnu, dtype=v.dtype)

        for i in range(nmu):
            mui = tuple(mu[i, :])
            j = lut[mui]
            res[j] = 1.0 / sqrt(multinomial(mui)) * v[i]

        return res


    def multiply_kronecker_power_v(self, D, A, v, i):
        r"""Multiply the :math:`i`-th Kronecker power of the matrix :math:`A`
        by a vector :math:`v` from the right.

        :param A: The matrix :math:`A`.
        :param v: The vector :math:`v`.
        :param i: The integer Kronecker power exponent.

        .. note:: The matrix has to be square.
        """
        assert i >= 0

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


    def multiply_T_v(self, coeffs, K, NU, MU, D, J):
        r"""Apply the transformation matrix :math:`T` to the coefficients :math:`c`.

        :param coeffs:
        :param K:
        :param NU:
        :param MU:
        :param D:
        :param J:
        """
        nDm = lambda D, m: int(binom(D + m - 1, m))

        s = [nDm(D, i) for i in range(J)]
        Ci = split(coeffs, cumsum(s))
        Cit = []

        for i, ci in enumerate(Ci):
            lut = self.build_lut(NU[i])
            cit = self.multiply_PiT_v(NU[i], MU[i], lut, ci)
            cit = self.multiply_kronecker_power_v(D, K, cit, i)
            cit = self.multiply_Pi_v(NU[i], MU[i], lut, cit)
            Cit.append(cit)

        return hstack(Cit).reshape(*coeffs.shape)


    def multiply_Tinv_v(self, coeffs, K, NU, MU, D, J):
        r"""Apply the transformation matrix :math:`T` to the coefficients :math:`d`.

        :param coeffs:
        :param K:
        :param NU:
        :param MU:
        :param D:
        :param J:
        """
        nDm = lambda D, m: int(binom(D + m - 1, m))

        s = [nDm(D, i) for i in range(J)]
        Ci = split(coeffs, cumsum(s))
        Cit = []

        Kinv = K.transpose().conjugate()

        for i, ci in enumerate(Ci):
            lut = self.build_lut(NU[i])
            cit = self.multiply_PiT_v(NU[i], MU[i], lut, ci)
            cit = self.multiply_kronecker_power_v(D, Kinv, cit, i)
            cit = self.multiply_Pi_v(NU[i], MU[i], lut, cit)
            Cit.append(cit)

        return hstack(Cit).reshape(*coeffs.shape)


    def _transform(self, HAWPfrom, LOP, HAWPtoConstructor):
        r"""The transformation logic. This function can transform
        in both directions depending on the linear operator argument.

        :param HAWPfrom: The wavepacket to transform.
        :param LOP: The linear transformation operator, either :math:`T` or :math:`T^{-1}`.
        :param HAWPtoConstructor: Constructor of the resulting wavepacket type.
        """
        if not HAWPfrom.get_number_components() == 1:
            raise NotImplementedError("Only scalar wavepackets are supported.")

        D = HAWPfrom.get_dimension()
        BS = HAWPfrom.get_basis_shapes(component=0)

        if not isinstance(BS, SimplexShape):
            raise ValueError("The wavepacket does not have a simplex basis shape.")

        J = BS.get_description()['K']

        NU = self.build_nu(D, J)
        MU = self.build_mu(D, J)

        # Compute small overlap matrix
        Pi = HAWPfrom.get_parameters()
        eps = HAWPfrom.get_eps()
        M = self.overlap(D, Pi, eps)

        # Permutation, adapter to basis index order
        cfrom = HAWPfrom.get_coefficients(component=0)
        cfromt = self.mapitinv(BS, MU, cfrom)
        ctot = LOP(cfromt, M, NU, MU, D, J)
        cto = self.mapit(BS, MU, ctot)

        # Set up a to wavepacket
        HAWPto = HAWPtoConstructor(D, 1, eps)
        HAWPto.set_parameters(Pi)
        HAWPto.set_basis_shapes(BS, component=0)
        HAWPto.set_coefficients(cto, component=0)

        return HAWPto


    def transform_phi_to_psi(self, HAWPphi):
        r"""Transform a old-kind wavepacket :math:`\Phi` into a new-kind wavepacket :math:`\Psi`.

        :param HAWPphi: The wavepacket :math:`\Phi` to transform.
        :return: A new wavepacket object for :math:`\Psi`.
        """
        return self._transform(HAWPphi, self.multiply_T_v, HagedornWavepacketNew)


    def transform_psi_to_phi(self, HAWPpsi):
        r"""Transform a new-kind wavepacket :math:`\Psi` into a old-kind wavepacket :math:`\Phi`.

        :param HAWPphi: The wavepacket :math:`\Psi` to transform.
        :return: A new wavepacket object for :math:`\Phi`.
        """
        return self._transform(HAWPpsi, self.multiply_Tinv_v, HagedornWavepacket)
