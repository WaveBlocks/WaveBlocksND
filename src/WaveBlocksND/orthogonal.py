"""The WaveBlocks Project

This file implements computation of Guass-Hermite
quadrature nodes and weights. For high order a novel
asymptotic linear time algorithm is used.

@author: R. Bourquin
@copyright: Copyright (C) 2014, 2015 R. Bourquin
@license: Modified BSD License
"""

from numpy import (sqrt, exp, sin, cos, arccos, pi, floor, ceil, around, arange,
                   array, ones_like, zeros_like, row_stack, column_stack, hstack,
                   flipud, floating, dot)
from scipy.special import h_roots as h_roots_original, airy


def compute_tauk(n, k, maxit=5):
    """Helper function for Tricomi initial guesses

    For details, see formula 3.1 in lemma 3.1 in the
    original paper.

    Parameters
    ----------
    n : int
        Quadrature order
    k : ndarray of type int
        Index of roots :math:`\tau_k` to compute
    maxit : int
        Number of Newton maxit performed, the default
        value of 5 is sufficient.

    Returns
    -------
    tauk : ndarray
        Roots of equation 3.1

    See Also
    --------
    initial_nodes_a
    h_roots_asy
    """
    a = n % 2 - 0.5
    c = (4.0*floor(n/2.0) - 4.0*k + 3.0)*pi / (4.0*floor(n/2.0) + 2.0*a + 2.0)
    f = lambda x: x - sin(x) - c
    df = lambda x: 1.0 - cos(x)
    xi = 0.5*pi
    for i in range(maxit):
        xi = xi - f(xi)/df(xi)
    return xi


def initial_nodes_a(n, k):
    """Tricomi initial guesses

    Computes an initial approximation to the square of the `k`-th
    (positive) root :math:`x_k` of the Hermite polynomial :math:`H_n`
    of order :math:`n`. The formula is the one from lemma 3.1 in the
    original paper. The guesses are accurate except in the region
    near :math:`\sqrt{2n + 1}`.

    Parameters
    ----------
    n : int
        Quadrature order
    k : ndarray of type int
        Index of roots to compute

    Returns
    -------
    xksq : ndarray
        Square of the approximate roots

    See Also
    --------
    initial_nodes
    h_roots_asy
    """
    tauk = compute_tauk(n, k)
    sigk = cos(0.5*tauk)**2
    a = n % 2 - 0.5
    nu = 4.0*floor(n/2.0) + 2.0*a + 2.0
    # Initial approximation of Hermite roots (square)
    xksq = nu*sigk - 1.0/(3.0*nu) * (5.0/(4.0*(1.0-sigk)**2) - 1.0/(1.0-sigk) - 0.25)
    return xksq


def compute_am(m):
    """Approximation of the roots of the Airy function

    Computes an asymptotic approximation to the root :math:`a_m`
    of the Airy function :math:`Ai(x)` for large values of :math:`m`.
    This formula is essentially exact for :math:`m > 10`.

    Parameters
    ----------
    m : ndarray of type int
        Index of the root :math:`a_m`

    Returns
    -------
    am : ndarray
        The `m`-th root :math:`a_m` of the function :math:`Ai(x)`

    See Also
    --------
    airy_root
    h_roots_asy
    """
    sm = 3*pi*(4*m-1) / 8.0
    coeffs = row_stack([1.0, 5.0/48.0, -5.0/36.0, 77125.0/82944.0,
                        -108056875.0/6967296.0, 162375596875.0/334430208.0])
    smp = column_stack([ones_like(sm), sm**(-2), sm**(-4), sm**(-6), sm**(-8), sm**(-10)])
    am = -sm**(2.0/3.0) * dot(smp, coeffs).reshape((-1,))
    return am


def airy_root(m):
    """Approximation of the roots of the Airy function

    Computes an approximation to the `m`-th root :math:`a_m` of the
    Airy function :math:`Ai(x)`. For :math:`m < 10` we use a lookup
    table and an asymptotic expansion otherwise.

    Parameters
    ----------
    m : ndarray of type int
        Index of the root :math:`a_m`

    Returns
    -------
    am : ndarray
        The `m`-th root :math:`a_m` of the function :math:`Ai(x)`

    See Also
    --------
    h_roots_asy
    """
    airyroots = array([
             -2.3381074104597670385,
             -4.0879494441309706166,
             -5.5205598280955510591,
             -6.7867080900717589988,
             -7.9441335871208531231,
             -9.0226508533409803802,
             -10.040174341558085931,
             -11.008524303733262893,
             -11.936015563236262517,
             -12.828776752865757200])
    # Note: This specialized code is much faster than 'ai_zeros'
    am = zeros_like(m, dtype=floating)
    I = m <= 10
    am[I] = airyroots[m[I]-1]
    J = m > 10
    am[J] = compute_am(m[J])
    return am


def initial_nodes_b(n, k):
    """Gatteschi initial guesses

    Computes an initial approximation to the square of the `k`-th
    (positive) root :math:`x_k` of the Hermite polynomial :math:`H_n`
    of order :math:`n`. The formula is the one from lemma 3.2 in the
    original paper. The guesses are accurate in the region just
    below :math:`\sqrt{2n + 1}`.

    Parameters
    ----------
    n : int
        Quadrature order
    k : ndarray of type int
        Index of roots to compute

    Returns
    -------
    xksq : ndarray
        Square of the approximate root

    See Also
    --------
    initial_nodes
    h_roots_asy
    """
    a = n % 2 - 0.5
    nu = 4.0*floor(n/2.0) + 2.0*a + 2.0
    # Airy roots by approximation
    ak = airy_root(k)
    # Initial approximation of Hermite roots (square)
    xksq = (nu +
            2.0**(2.0/3.0) * ak * nu**(1.0/3.0) +
            1.0/5.0 * 2.0**(4.0/3.0) * ak**2 * nu**(-1.0/3.0) +
            (9.0/140.0 - 12.0/175.0 * ak**3) * nu**(-1.0) +
            (16.0/1575.0 * ak + 92.0/7875.0 * ak**4) * 2.0**(2.0/3.0) * nu**(-5.0/3.0) -
            (15152.0/3031875.0 * ak**5 + 1088.0/121275.0 * ak**2) * 2.0**(1.0/3.0) * nu**(-7.0/3.0))
    return xksq


def initial_nodes(n):
    """Initial guesses for the Hermite roots

    Computes an initial approximation to the non-negative
    roots :math:`x_k` of the Hermite polynomial :math:`H_n`
    of order :math:`n`. The Tricomi and Gatteschi initial
    guesses are used in the region where they are accurate.

    Parameters
    ----------
    n : int
        Quadrature order

    Returns
    -------
    xk : ndarray
        Approximate roots

    See Also
    --------
    h_roots_asy
    """
    # Turnover point
    # linear polynomial fit to error of 10, 25, 40, ..., 1000 point rules
    fit = 0.49082003*n - 4.37859653
    turnover = around(fit).astype(int)
    # Compute all approximations
    ia = arange(1, int(floor(n*0.5)+1))
    ib = flipud(arange(1, int(1+n-ceil(n*0.5))))
    xasq = initial_nodes_a(n, ia[:turnover+1])
    xbsq = initial_nodes_b(n, ib[turnover+1:])
    # Combine
    iv = sqrt(hstack([xasq, xbsq]))
    # Central node is always zero
    if n % 2 == 1:
        iv = hstack([0.0, iv])
    return iv


def pbcf(n, theta):
    """Asymptotic series expansion of parabolic cylinder function

    The implementation is based on sections 3.2 and 3.3 from the
    original paper. Compared to the published version this code
    adds one more term to the asymptotic series. The detailed
    formulas can be found at [parabolic-asymptotics]_. The evaluation
    is done in a transformed variable :math:`\theta := \arccos(t)`
    where :math:`t := x / \mu` and :math:`\mu := \sqrt{2n + 1}`.

    Parameters
    ----------
    n : int
        Quadrature order
    theta : ndarray
        Transformed position variable

    Returns
    -------
    U : ndarray
        Value of the parabolic cylinder function :math:`U(a, \theta)`.
    Ud : ndarray
        Value of the derivative :math:`U^{\prime}(a, \theta)` of
        the parabolic cylinder function.

    See Also
    --------
    h_roots_asy

    References
    ----------
    .. [parabolic-asymptotics]
       http://dlmf.nist.gov/12.10#vii
    """
    mu = sqrt(2.0*n + 1.0)
    st = sin(theta)
    ct = cos(theta)
    # http://dlmf.nist.gov/12.10#vii
    mu = 2.0*n + 1.0
    # http://dlmf.nist.gov/12.10#E23
    eta = 0.5*theta - 0.5*st*ct
    # http://dlmf.nist.gov/12.10#E39
    zeta = -(3.0*eta/2.0) ** (2.0/3.0)
    # http://dlmf.nist.gov/12.10#E40
    phi = (-zeta / st**2) ** (0.25)
    # Coefficients
    # http://dlmf.nist.gov/12.10#E43
    a0 =  1.0
    a1 =  0.10416666666666666667
    a2 =  0.08355034722222222222
    a3 =  0.12822657455632716049
    a4 =  0.29184902646414046425
    a5 =  0.88162726744375765242
    b0 =  1.0
    b1 = -0.14583333333333333333
    b2 = -0.09874131944444444444
    b3 = -0.14331205391589506173
    b4 = -0.31722720267841354810
    b5 = -0.94242914795712024914
    a = (a0, a1, a2, a3, a4, a5)
    b = (b0, b1, b2, b3, b4, b5)
    # Polynomials
    # http://dlmf.nist.gov/12.10#E9
    # http://dlmf.nist.gov/12.10#E10
    u0 = 1.0
    ctp = ct ** arange(16).reshape((-1,1))
    u1 = (       1.0*ctp[3,:] -          6.0*ct) / 24.0
    u2 = (      -9.0*ctp[4,:] +        249.0*ctp[2,:] +         145.0) / 1152.0
    u3 = (   -4042.0*ctp[9,:] +      18189.0*ctp[7,:] -       28287.0*ctp[5,:] -      151995.0*ctp[3,:] -     259290.0*ct) / 414720.0
    u4 = (   72756.0*ctp[10,:] -    321339.0*ctp[8,:] -      154982.0*ctp[6,:] +    50938215.0*ctp[4,:] +  122602962.0*ctp[2,:] + 12773113.0) / 39813120.0
    u5 = (82393456.0*ctp[15,:] - 617950920.0*ctp[13,:] + 1994971575.0*ctp[11,:] - 3630137104.0*ctp[9,:] + 4433574213.0*ctp[7,:]
          - 37370295816.0*ctp[5,:] - 119582875013.0*ctp[3,:] - 34009066266.0*ct) / 6688604160.0
    v0 = 1.0
    v1 = (       1.0*ctp[3,:] +           6.0*ct) / 24.0
    v2 = (      15.0*ctp[4,:] -         327.0*ctp[2,:] -         143.0) / 1152.0
    v3 = (   -4042.0*ctp[9,:] +       18189.0*ctp[7,:] -       36387.0*ctp[5,:] +      238425.0*ctp[3,:] +     259290.0*ct) / 414720.0
    v4 = ( -121260.0*ctp[10,:] +     551733.0*ctp[8,:] -      151958.0*ctp[6,:] -    57484425.0*ctp[4,:] -  132752238.0*ctp[2,:] - 12118727) / 39813120.0
    v5 = (82393456.0*ctp[15,:] - 617950920.0*ctp[13,:] + 2025529095.0*ctp[11,:] - 3750839308.0*ctp[9,:] + 3832454253.0*ctp[7,:]
          + 35213253348.0*ctp[5,:] + 130919230435.0*ctp[3,:] + 34009066266*ct) / 6688604160.0
    u = (u0, u1, u2, u3, u4, u5)
    v = (v0, v1, v2, v3, v4, v5)
    # Airy Evaluation (Bi and Bip unused)
    Ai, Aip, Bi, Bip = airy(mu**(4.0/6.0) * zeta)
    # Prefactor for U
    P = 2.0*sqrt(pi) * mu**(1.0/6.0) * phi
    # Terms for U
    # http://dlmf.nist.gov/12.10#E42
    A0 =   b[0]*u[0]
    A1 =  (b[2]*u[0] + phi**6*b[1]*u[1] + phi**12*b[0]*u[2]) / zeta**3
    A2 =  (b[4]*u[0] + phi**6*b[3]*u[1] + phi**12*b[2]*u[2] + phi**18*b[1]*u[3] + phi**24*b[0]*u[4]) / zeta**6
    B0 = -(a[1]*u[0] + phi**6*a[0]*u[1]) / zeta**2
    B1 = -(a[3]*u[0] + phi**6*a[2]*u[1] + phi**12*a[1]*u[2] + phi**18*a[0]*u[3]) / zeta**5
    B2 = -(a[5]*u[0] + phi**6*a[4]*u[1] + phi**12*a[3]*u[2] + phi**18*a[2]*u[3] + phi**24*a[1]*u[4] + phi**30*a[0]*u[5]) / zeta**8
    # U
    # http://dlmf.nist.gov/12.10#E35
    U = P * (Ai  * (A0 + A1/mu**2.0 + A2/mu**4.0) +
             Aip * (B0 + B1/mu**2.0 + B2/mu**4.0) / mu**(8.0/6.0))
    # Prefactor for derivative of U
    Pd = sqrt(2.0*pi) * mu**(2.0/6.0) / phi
    # Terms for derivative of U
    # http://dlmf.nist.gov/12.10#E46
    C0 = -(b[1]*v[0] + phi**6*b[0]*v[1]) / zeta
    C1 = -(b[3]*v[0] + phi**6*b[2]*v[1] + phi**12*b[1]*v[2] + phi**18*b[0]*v[3]) / zeta**4
    C2 = -(b[5]*v[0] + phi**6*b[4]*v[1] + phi**12*b[3]*v[2] + phi**18*b[2]*v[3] + phi**24*b[1]*v[4] + phi**30*b[0]*v[5]) / zeta**7
    D0 =   a[0]*v[0]
    D1 =  (a[2]*v[0] + phi**6*a[1]*v[1] + phi**12*a[0]*v[2]) / zeta**3
    D2 =  (a[4]*v[0] + phi**6*a[3]*v[1] + phi**12*a[2]*v[2] + phi**18*a[1]*v[3] + phi**24*a[0]*v[4]) / zeta**6
    # Derivative of U
    # http://dlmf.nist.gov/12.10#E36
    Ud = Pd * (Ai  * (C0 + C1/mu**2.0 + C2/mu**4.0) / mu**(4.0/6.0) +
               Aip * (D0 + D1/mu**2.0 + D2/mu**4.0) )
    return U, Ud


def newton(n, x_initial, maxit=5):
    """Newton iteration for polishing the asymptotic approximation
    to the zeros of the Hermite polynomials.

    Parameters
    ----------
    n : int
        Quadrature order
    x_initial : ndarray
        Initial guesses for the roots
    maxit : int
        Maximal number of Newton iterations.
        The default 5 is sufficient, usually
        only one or two steps are needed.

    Returns
    -------
    nodes : ndarray
        Quadrature nodes
    weights : ndarray
        Quadrature weights

    See Also
    --------
    h_roots_asy
    """
    # Variable transformation
    mu = sqrt(2.0*n + 1.0)
    t = x_initial / mu
    theta = arccos(t)
    # Newton iteration
    for i in range(maxit):
        u, ud = pbcf(n, theta)
        dtheta = u / (sqrt(2.0) * mu * sin(theta) * ud)
        theta = theta + dtheta
        if max(abs(dtheta)) < 1e-14:
            break
    # Undo variable transformation
    x = mu * cos(theta)
    # Central node is always zero
    if n % 2 == 1:
        x[0] = 0.0
    # Compute weights
    w = 1.0 / (2.0*ud**2)
    return x, w


def h_roots_asy(n):
    """Gauss-Hermite (physicst's) quadrature for large n

    Computes the sample points and weights for Gauss-Hermite quadrature.
    The sample points are the roots of the `n`th degree Hermite polynomial,
    :math:`H_n(x)`.  These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-inf, inf]` with weight function :math:`f(x) = e^{-x^2}`.
    This method relies on asymptotic expansions which work best for n > 150.
    The algorithm has linear runtime making computation for very large n
    feasible.

    Parameters
    ----------
    n : int
        quadrature order

    Returns
    -------
    nodes : ndarray
        Quadrature nodes
    weights : ndarray
        Quadrature weights

    See Also
    --------
    h_roots

    References
    ----------
    .. [townsend.trogdon.olver-2014]
       Townsend, A. and Trogdon, T. and Olver, S. (2014)
       *Fast computation of Gauss quadrature nodes and
       weights on the whole real line*. ArXiv 1410.5286.
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    iv = initial_nodes(n)
    nodes, weights = newton(n, iv)
    # Combine with negative parts
    if n % 2 == 0:
        nodes = hstack([-flipud(nodes), nodes])
        weights = hstack([flipud(weights), weights])
    else:
        nodes = hstack([-flipud(nodes[1:]), nodes])
        weights = hstack([flipud(weights[1:]), weights])
    # Scale weights
    weights *= exp(-nodes**2)
    weights *= sqrt(pi) / sum(weights)
    return nodes, weights


def h_roots(n):
    """Gauss-Hermite (physicst's) quadrature

    Computes the sample points and weights for Gauss-Hermite quadrature.
    The sample points are the roots of the `n`th degree Hermite polynomial,
    :math:`H_n(x)`. These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-inf, inf]` with weight function :math:`f(x) = e^{-x^2}`.

    Parameters
    ----------
    n : int
        quadrature order

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights

    Notes
    -----
    For small n up to 150 a modified version of the Golub-Welsch
    algorithm is used. Nodes are computed from the eigenvalue
    problem and improved by one step of a Newton iteration.
    The weights are computed from the well-known analytical formula.

    For n larger than 150 an optimal asymptotic algorithm is applied
    which computes nodes and weights in a numerically stable manner.
    The algorithm has linear runtime making computation for very
    large n (several thousand or more) feasible.

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.hermite.hermgauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    if n <= 150:
        return h_roots_original(n)
    else:
        return h_roots_asy(n)


def hermite_recursion(n, x):
    r"""Evaluate the Hermite function :math:`h_n(x)` of order :math:`n`
    recursively on the given points :math:`x` by the three term recursion.

    :param n: The order :math:`n` of the Hermite function :math:`h_n(x)` to evaluate.
    :param x: The points at which the Hermite function :math:`h_n(x)` is evaluated.
    :return: An array :math:`Hn` containing the values of :math:`h_n(x)`.
    """
    Hnm1 = zeros_like(x)
    Hn   = zeros_like(x)
    Hnp1 = zeros_like(x)

    Hn = pi**(-0.25) * exp(-0.5*x**2)
    if n >= 1:
        Hnm1 = Hn
        Hn = sqrt(2.0) * x * Hnm1
        for k in xrange(2, n+1):
            Hnp1 = sqrt(2.0/k) * x * Hn - sqrt((k-1.0)/k) * Hnm1
            Hnm1 = Hn
            Hn = Hnp1

    return Hn


def psi_roots_asy(n):
    """Gauss-Hermite (physicst's) quadrature for large n

    Computes the sample points and weights for transformed Gauss-Hermite
    quadrature. The sample points are the roots of the `n`th degree Hermite
    polynomial, :math:`H_n(x)`. The weights are transformed such that
    they do not include the exponential factor :math:`\exp(-x^2)`.

    Parameters
    ----------
    n : int
        quadrature order

    Returns
    -------
    nodes : ndarray
        Quadrature nodes
    weights : ndarray
        Quadrature weights

    See Also
    --------
    psi_roots
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    iv = initial_nodes(n)
    nodes, weights = newton(n, iv)
    # Combine with negative parts
    if n % 2 == 0:
        nodes = hstack([-flipud(nodes), nodes])
        weights = hstack([flipud(weights), weights])
    else:
        nodes = hstack([-flipud(nodes[1:]), nodes])
        weights = hstack([flipud(weights[1:]), weights])
    # Scale weights
    C = sqrt(2)*sqrt(pi) / sum(exp(-0.5*nodes**2) * weights)
    weights *= C
    return nodes, weights


def psi_roots(n):
    """Gauss-Hermite quadrature for functions including the exponential factor.

    Computes the sample points and weights for transformed Gauss-Hermite
    quadrature. The sample points are the roots of the `n`th degree Hermite
    polynomial, :math:`H_n(x)`. The weights are transformed such that
    they do not include the exponential factor :math:`\exp(-x^2)`.

    Parameters
    ----------
    n : int
        quadrature order

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights

    Notes
    -----
    For small n up to 150 a modified version of the Golub-Welsch
    algorithm is used. Nodes are computed from the eigenvalue
    problem and improved by one step of a Newton iteration.
    The weights are computed from the well-known analytical formula.

    For n larger than 150 an optimal asymptotic algorithm is applied
    which computes nodes and weights in a numerically stable manner.
    The algorithm has linear runtime making computation for very
    large n (several thousand or more) feasible.

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.hermite.hermgauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    if n <= 150:
        nodes, weights =  h_roots_original(n)
        # Hermite function recursion
        h = hermite_recursion(n-1, nodes)
        weights = 1.0/(h**2 * n)
        return nodes, weights
    else:
        return psi_roots_asy(n)
