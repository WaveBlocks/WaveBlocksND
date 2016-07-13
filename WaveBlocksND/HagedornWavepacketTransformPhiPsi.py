"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from collections import defaultdict

from numpy import (complexfloating, array, hstack, vstack, split, integer,
                   conjugate, zeros, identity, zeros_like, cumsum)
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


    def mapit(BS, MU, c):
        # TODO: Will resolve in identity
        cm = zeros_like(c)
        K = BS.get_basis_size()
        mu = vstack(MU)
        for k in range(K):
            cm[BS[tuple(mu[k,:])]] = c[k]
        return cm


    def mapitinv(BS, MU, c):
        # TODO: Will resolve in identity
        cm = zeros_like(c)
        K = BS.get_basis_size()
        mu = vstack(MU)
        for k in range(K):
            cm[k] = c[BS[tuple(mu[k,:])]]
        return cm


    def overlap(self, D, Pi, eps):
        r"""
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
        r"""
        Build the lookup table for simultaneuos multi-acces to arrays by labellists.
        """
        lut = defaultdict(list)

        for i, k in enumerate(nu):
            lut[tuple(k)].append(i)

        # TODO: Maybe convert lists to arrays?
        return lut


    def multiply_Pi_v(self, nu, mu, lut, v):
        r"""
        Multiply the matrix :math:`P_i` by a vector :math:`v` from the right.

        :param nu: The list :math:`\nu`.
        :param mu: The list :math:`\mu`.
        :param lut: The lookup table.
        :param v: The vector :math:`v \in \mathcal{X}_i`.
        """
        nmu = mu.shape[0]
        res = zeros((nmu), dtype=v.dtype)

        for i in range(nmu):
            mui = tuple(mu[i, :])
            j = lut[mui][0]
            res[i] = sqrt(multinomial(mui)) * v[j]

        return res


    def multiply_PiT_v(self, nu, mu, lut, v):
        r"""
        Multiply the matrix :math:`P_i^{T}` by a vector :math:`v` from the right.

        :param nu: The list :math:`\nu`.
        :param mu: The list :math:`\mu`.
        :param lut: The lookup table.
        :param v: The vector :math:`v \in \mathbb{C}^{n(D,i)}`.
        """
        nmu = mu.shape[0]
        nnu = nu.shape[0]
        res = zeros((nnu), dtype=v.dtype)

        for i in range(nmu):
            mui = tuple(mu[i, :])
            j = lut[mui]
            res[j] = 1.0 / sqrt(multinomial(mui)) * v[i]

        return res


    def multiply_kronecker_power_v(self, A, v, i):
        r"""
        Multiply the :math:`i`-th Kronecker power of the matrix :math:`A` by a
        vector :math:`v` from the right.

        :param A: The matrix :math:`A`.
        :param v: The vector :math:`v`.
        :param i: The integer Kronecker power exponent.

        .. note:: The matrix has to be square.
        """
        assert i >= 0

        select = lambda data, base, size, step: data.reshape(-1, step)[:, base*size:(base+1)*size]

        if i == 0:
            return v

        assert A.shape[0] == A.shape[1]
        D = A.shape[0]

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
        r"""
        Apply the transformation matrix :math:`T` to the coefficients :math:`c`.

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
            cit = self.multiply_kronecker_power_v(K, cit, i)
            cit = self.multiply_Pi_v(NU[i], MU[i], lut, cit)
            Cit.append(cit)

        return hstack(Cit).reshape(*coeffs.shape)


    def multiply_Tinv_v(self, coeffs, K, NU, MU, D, J):
        r"""
        """
        nDm = lambda D, m: int(binom(D + m - 1, m))

        raise NotImplementedError()


    def transform_phi_psi(self, HAWPold):
        r"""
        # Old -> New
        # Muliply by T
        """
        D = HAWPold.get_dimension()

        # TODO: We handle scalar wavepackets only
        BS = HAWPold.get_basis_shapes(component=0)

        if not isinstance(BS, SimplexShape):
            raise ValueError("The wavepacket does not have a simplex basis shape.")

        # J = sum(BS.find_largest_index())
        J = BS.get_description()['K']

        NU = self.build_nu(D, J)
        MU = self.build_mu(D, J)

        # Compute small overlap matrix
        Pi = HAWPold.get_parameters()
        eps = HAWPold.get_eps()

        M = self.overlap(D, Pi, eps)

        # Permutation, adapter to basis index order
        cold = HAWPold.get_coefficients(component=0)
        coldt = self.mapitinv(BS, MU, cold)
        cnewt = self.multiply_T_v(coldt, M, NU, MU, D, J)
        cnew = self.mapit(BS, MU, cnewt)

        # Set up a new wavepacket
        HAWPnew = HagedornWavepacketNew(D, 1, eps)
        HAWPnew.set_parameters(Pi)
        HAWPnew.set_basis_shapes([BS])
        HAWPnew.set_coefficients(cnew)

        return HAWPnew


    def transform_psi_phi(self, HAWPnew):
        r"""
        # New -> Old
        # Multiply by T^{-1}
        """
        D = HAWPnew.get_dimension()

        # TODO: We handle scalar wavepackets only
        BS = HAWPnew.get_basis_shapes(component=0)

        if not isinstance(BS, SimplexShape):
            raise ValueError("The wavepacket does not have a simplex basis shape.")

        # J = sum(BS.find_largest_index())
        J = BS.get_description()['K']

        NU = self.build_nu(D, J)
        MU = self.build_mu(D, J)

        # Compute small overlap matrix
        Pi = HAWPnew.get_parameters()
        eps = HAWPnew.get_eps()

        M = self.overlap(D, Pi, eps)

        # Permutation, adapter to basis index order
        cnew = HAWPnew.get_coefficients(component=0)
        cnewt = self.mapitinv(BS, MU, cnew)
        coldt = self.multiply_Tinv_v(cnewt, M, NU, MU, D, J)
        cold = self.mapit(BS, MU, coldt)

        # Set up a new wavepacket
        HAWPold = HagedornWavepacketNew(D, 1, eps)
        HAWPold.set_parameters(Pi)
        HAWPold.set_basis_shapes([BS])
        HAWPold.set_coefficients(cold)

        return HAWPold
