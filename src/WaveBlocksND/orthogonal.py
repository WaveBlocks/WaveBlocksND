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
from scipy.special.orthogonal import _h_gen_roots_and_weights
from scipy.special.orthogonal import cephes


def _h_gen_roots_and_weights_large(n, mu, factor, func):
    r"""Compute the roots and weights for Gaussian-Hermite quadrature.
    Internal function.
    """
    if n < 1:
        raise ValueError("n must be positive.")

    bn = np.sqrt(np.arange(1, n, dtype=np.float64)/factor)
    c = np.diag(bn, -1)
    x, ev = linalg.eigh(c)
    w = ev[0,:]**2
    # symmetrize
    w = (w + w[::-1])/2
    x = (x - x[::-1])/2
    # scale w correctly
    w *= np.sqrt(2.0*np.pi/factor)
    if mu:
        return [x, w, mu]
    else:
        return x, w


def h_roots(n, mu=0):
    r"""[x,w] = h_roots(n)

    Returns the roots (x) of the nth order Hermite polynomial,
    H_n(x), and weights (w) to use in Gaussian Quadrature over
    [-inf,inf] with weighting function exp(-x**2).
    """
    n_max = 200
    if n <= n_max:
        return _h_gen_roots_and_weights(n, mu, 2.0, cephes.eval_hermite)
    else:
        return _h_gen_roots_and_weights_large(n, mu, 2.0, cephes.eval_hermite)
