"""The WaveBlocks Project

This file contains several different algorithms to compute the
matrix exponential. Currently we have an exponential based on
Pade approximations and an Arnoldi iteration method.

@author: R. Bourquin
@copyright: Copyright (C) 2007 V. Gradinaru
@copyright: Copyright (C) 2010, 2011, 2012, 2015 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, dot, complexfloating, conjugate
from scipy.linalg import norm, expm


def matrix_exp_pade(A, v, factor):
    r"""Compute the solution of :math:`v' = A v` with a full
    matrix exponential via Pade approximation.

    :param A: The matrix :math:`A` of shape :math:`N \times N`.
    :param v: The vector :math:`v` of length :math:`N`.
    :param factor: An additional scalar factor :math:`\alpha`.
    :return: The (approximate) value of :math:`\exp\left(\alpha A\right) v`
    """
    return dot(expm(factor * A), v)


def arnoldi(A, v0, k):
    r"""Arnoldi algorithm to compute the Krylov approximation :math:`H` of a matrix :math:`A`.

    :param A: The matrix :math:`A` of shape :math:`N \times N` to approximate.
    :param v0: The initial vector :math:`v_0` of length :math:`N`.
    :param k: The number :math:`k` of Krylov steps performed.
    :return: A tuple :math:`(V, H)` where :math:`V` is the large matrix of shape
             :math:`N \times (k+1)` containing the orthogonal vectors and :math:`H` is the
             small matrix of shape :math:`(k+1) \times k` containing the Krylov approximation
             of :math:`A`.
    """
    r, c = A.shape
    V = zeros((r, k + 1), dtype=complexfloating)
    H = zeros((k + 1, k), dtype=complexfloating)

    V[:, 0] = v0.reshape(-1) / norm(v0)

    for i in range(1, k + 1):
        vi = dot(A, V[:, i - 1])
        for j in range(i):
            H[j, i - 1] = dot(conjugate(V[:, j]), vi)
            vi -= H[j, i - 1] * V[:, j]
        H[i, i - 1] = norm(vi)
        V[:, i] = vi / H[i, i - 1]

    return V, H


def matrix_exp_arnoldi(A, v, factor, k):
    r"""Compute the solution of :math:`v' = A v` via :math:`k`
    steps of a the Arnoldi krylov method.

    :param A: The matrix :math:`A` of shape :math:`N \times N`.
    :param v: The vector :math:`v` of length :math:`N`.
    :param factor: An additional scalar factor :math:`\alpha`.
    :param k: The number :math:`k` of Krylov steps performed.
    :return: The (approximate) value of :math:`\exp\left(\alpha A\right) v`.
    """
    V, H = arnoldi(A, v, min(min(A.shape), k))
    eH = expm(factor * H[:-1, :])
    r = norm(v) * dot(V[:, :-1], eH[:, 0])
    return r.reshape(v.shape)
