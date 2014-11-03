"""The WaveBlocks Project

Combinatorics functions for enumerating various objects.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, integer


def partitions(D, K):
    r"""Enumerate integer partitions in anti-lexicographic
    order for integers up to (and including) some limit K.
    All partitions have exactly D parts, some may be zero.

    :param D: Dimension
    :param K: Limit
    """
    P = zeros((D,), dtype=integer)

    yield P.copy()
    while P.sum() <= K:
        p0 = P[0]
        for i in xrange(1, D):
            p0 += P[i]
            if P[0] <= P[i] + 1:
                P[i] = 0
            else:
                P[0] = p0 - i * (P[i] + 1)
                P[1:(i+1)] = P[i] + 1
                yield P.copy()
                break
        else:
            P[0] = p0 + 1
            if P.sum() <= K:
                assert P.sum() <= K
                yield P.copy()


def lattice_points(D, N):
    r"""Enumerate all lattice points of an integer lattice
    :math:`\Lambda \subset \mathbb{N}^D` in :math:`D` dimensions
    having :math:`l_1` norm smaller or equal to :math:`N`.

    :param D: The dimension :math:`D` of the lattice.
    :param N: The maximal :math:`l_1` norm of the lattice points.
    """
    for n in xrange(N+1):
        k = zeros(D, dtype=integer)
        k[0] = n
        yield tuple(k)
        c = 1
        while k[D-1] < n:
            if c == D:
                for i in xrange(c-1, 0, -1):
                    c = i
                    if not k[i-1] == 0:
                        break
            k[c-1] = k[c-1] - 1
            c += 1
            k[c-1] = n - sum(k[0:c-1])
            if c < D:
                k[c:D] = zeros(D-c, dtype=integer)
            yield tuple(k)


def permutations(P):
    r"""Enumerate all permutations in anti-lexicographical
    order follwing the given permutation `P`.

    :param P: A permutation
    :type P: An `ndarray` of shape `(D,)`
    """
    D = P.size
    P = P.copy()
    yield P.copy()
    while True:
        for i in xrange(1, D):
            pi = P[i]
            if P[i-1] > pi:
                I = i
                if i > 1:
                    J = I
                    for j in xrange(I//2):
                        pj = P[j]
                        if pj <= pi:
                            I = I - 1
                        P[j] = P[i - j-1]
                        P[i - j-1] = pj
                        if P[j] > pi:
                            J = j+1
                    if P[I-1] <= pi:
                        I = J
                P[i] = P[I-1]
                P[I-1] = pi
                yield P.copy()
                break
        else:
            raise StopIteration
