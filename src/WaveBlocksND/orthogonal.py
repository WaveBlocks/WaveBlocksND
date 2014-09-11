"""The WaveBlocks Project

This file contains a patch for the 'h_roots' function from scipy.

The reason is that the scipy function 'h_roots' performes
a Newton step after computing the nodes via Golub-Welsch [1].
This fails for large enough n because a Hermite polynomial
of high order is evaluated.

This patch fixes the problem by not doing the Newton step
for large n. The resulting nodes and weights are not fully
correct but tiny anyway. In our usecase it does not matter.
Important is that we can avoid nan values appearing.

[1] BUG: special: improve h_roots accuracy at high order
    Commit: fadbd9ff177741aba57bde6a62cb909e0d4c5c36
    https://github.com/scipy/scipy/issues/1991

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

import numpy as np
from scipy import linalg
from scipy.special.orthogonal import cephes


def _gen_roots_and_weights(n, mu0, an_func, bn_func, f, df, weights_formula, symmetrize, mu):
    """[x,w] = gen_roots_and_weights(n,an_func,sqrt_bn_func,mu)

    Returns the roots (x) of an nth order orthogonal polynomial,
    and weights (w) to use in appropriate Gaussian quadrature with that
    orthogonal polynomial.

    The polynomials have the recurrence relation
          P_n+1(x) = (x - A_n) P_n(x) - B_n P_n-1(x)

    an_func(n)          should return A_n
    sqrt_bn_func(n)     should return sqrt(B_n)
    mu ( = h_0 )        is the integral of the weight over the orthogonal
                        interval
    """
    k = np.arange(n, dtype='d')
    c = np.zeros((2, n))
    c[0,1:] = bn_func(k[1:])
    c[1,:] = an_func(k)

    if weights_formula:
        x = linalg.eigvals_banded(c, overwrite_a_band=True)
        # improve roots by one application of Newton's method
        y = f(n, x)
        dy = df(n, x)
        x -= y/dy
        fm = f(n-1, x)
        fm /= np.abs(fm).max()
        dy /= np.abs(dy).max()
        w = 1.0 / (fm * dy)
    else:
        # Computing weights via explicit formula can fail for large n
        # because of the division by a polynomial of high degree.
        x, ev = linalg.eig_banded(c, overwrite_a_band=True)
        w = ev[0,:]**2

    if symmetrize:
        w = (w + w[::-1]) / 2
        x = (x - x[::-1]) / 2

    # scale w correctly
    w *= mu0 / w.sum()

    if mu:
        return x, w, mu0
    else:
        return x, w


def h_roots(n, mu=0):
    r"""[x,w] = h_roots(n)

    Returns the roots (x) of the nth order Hermite polynomial,
    H_n(x), and weights (w) to use in Gaussian Quadrature over
    [-inf,inf] with weighting function exp(-x**2).
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    mu0 = np.sqrt(np.pi)
    an_func = lambda k: 0.0*k
    bn_func = lambda k: np.sqrt(k/2.0)
    f = cephes.eval_hermite
    df = lambda n, x: 2.0 * n * cephes.eval_hermite(n-1, x)

    n_max = 200
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, n <= n_max, True, mu)
