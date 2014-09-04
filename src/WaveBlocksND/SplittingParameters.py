"""The WaveBlocks Project

This file contains data to build several closely
related splitting methods.

@author: V. Gradinaru, R. Bourquin
@copyright: Copyright (C) 2011, 2012, 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, flipud, diff

__all__ = ["SplittingParameters"]


class SplittingParameters(object):


    def build(self, method):
        u"""
        :param method: A string specifying the method for time integration.
        :return: Two arrays :math:`a` and :math:`b`.

        ====== ===== =========== =========
        Method Order Authors     Reference
        ====== ===== =========== =========
        LT     1     Lie/Trotter [1]_, [3]_ page 42, equation 5.2
        S2     2     Strang      [2]_, [3]_ page 42, equation 5.3
        SS     2     Strang      [2]_, [3]_ page 42, equation 5.3
        PRKS6  4     Blanes/Moan [4]_ page 318, table 2, 'S6'
        BM42   4     Blanes/Moan [4]_ page 318, table 3, 'SRKNb6'
        Y4     4     Yoshida     [5]_, [3]_ page 40, equation 4.4
        Y61    6     Yoshida     [5]_, [3]_ page 144, equation 3.11
        BM63   6     Blanes/Moan [4]_ page 318, table 3, 'SRKNa14'
        KL6    6     Kahan/Li    [6]_, [3]_ page 144, equation 3.12
        KL8    8     Kahan/Li    [6]_, [3]_ page 145, equation 3.14
        KL10   10    Kahan/Li    [6]_, [3]_ page 146, equation 3.15
        ====== ===== =========== =========

        .. [1] H.F. Trotter, "On the product of semi-groups of operators",
               Proc. Am. Math. Soc.1O (1959) 545-551.

        .. [2] G. Strang, "On the construction and comparison of difference schemes",
               SIAM J. Numer. Anal. 5 (1968) 506-517.

        .. [3] E. Hairer, C. Lubich, and G. Wanner, "Geometric Numerical Integration -
               Structure-Preserving Algorithms for Ordinary Differential Equations",
               Springer-Verlag, New York, 2002 (Corrected second printing 2004).

        .. [4] S. Blanes and P.C. Moan, "Practical Symplectic Partitioned
               Runge-Kutta and Runge-Kutta-Nystrom Methods", J. Computational and
               Applied Mathematics, Volume 142, Issue 2, (2002) 313-330.

        .. [5] H. Yoshida, "Construction of higher order symplectic integrators",
               Phys. Lett. A 150 (1990) 262-268.

        .. [6] W. Kahan and  R.-c. Li, "Composition constants for raising the orders
               of unconventional schemes for ordinary differential equations",
               Math. Comput. 66 (1997) 1089-1099.
        """
        if method == "LT":
            s = 1
            a = zeros(s)
            b = zeros(s)
            a[0] = 1.0
            b[0] = 1.0
        elif method == "S2":
            s = 2
            a = zeros(s)
            b = zeros(s)
            a[1] = 1.0
            b[0] = 0.5
            b[1] = 0.5
        elif method == "SS":
            s = 2
            a = zeros(s)
            b = zeros(s)
            a[0] = 0.5
            a[1] = 0.5
            b[0] = 1.0
        elif method =="PRKS6":
            s = 7
            a = zeros(s)
            b = zeros(s)
            a[1] =  0.209515106613362
            a[2] = -0.143851773179818
            a[3] =  0.5 - a[:3].sum()
            a[4:] = flipud(a[1:4])
            b[0] =  0.0792036964311957
            b[1] =  0.353172906049774
            b[2] = -0.0420650803577195
            b[3] =  1.0 - 2.0 * b[:3].sum()
            b[4:] = flipud(b[:3])
        elif method == "BM42":
            s = 7
            a = zeros(s)
            b = zeros(s)
            a[1] = 0.245298957184271
            a[2] = 0.604872665711080
            a[3] = 0.5 - a[:3].sum()
            a[4:] = flipud(a[1:4])
            b[0] =  0.0829844064174052
            b[1] =  0.396309801498368
            b[2] = -0.0390563049223486
            b[3] =  1.0 - 2.0 * b[:3].sum()
            b[4:] = flipud(b[:3])
        elif method == "Y4":
            s = 4
            a = zeros(s)
            b = zeros(s)
            pp = 3.0
            theta = 1.0/(2.0 - 2**(1.0/pp))
            vi = -2**(1.0/pp) * theta
            a[0] = 0.0
            a[1] = theta
            a[2] = vi
            a[3] = theta
            b[0] = 0.5 * theta
            b[1] = 0.5 * (vi + theta)
            b[2:] = flipud(b[:2])
        elif method == "Y61":
            s = 8
            a = zeros(s)
            b = zeros(s)
            a[1] =  0.78451361047755726381949763
            a[2] =  0.23557321335935813368479318
            a[3] = -1.17767998417887100694641568
            a[4] =  1.0 - 2.0 * a[1:4].sum()
            a[5:] = flipud(a[1:4])
            b[0] = 0.5 * a[1]
            b[1] = 0.5 * a[1:3].sum()
            b[2] = 0.5 * a[2:4].sum()
            b[3] = 0.5 * (1.0 - 4.0*b[1] - a[3])
            b[4:] = flipud(b[0:4])
        elif method == "BM63":
            s = 15
            a = zeros(s)
            b = zeros(s)
            a[1] =  0.09171915262446165
            a[2] =  0.183983170005006
            a[3] = -0.05653436583288827
            a[4] =  0.004914688774712854
            a[5] =  0.143761127168358
            a[6] =  0.328567693746804
            a[7] =  0.5 - a[:7].sum()
            a[8:] = flipud(a[1:8])
            b[0] =  0.0378593198406116
            b[1] =  0.102635633102435
            b[2] = -0.0258678882665587
            b[3] =  0.314241403071447
            b[4] = -0.130144459517415
            b[5] =  0.106417700369543
            b[6] = -0.00879424312851058
            b[7] =  1.0 - 2.0 * b[:7].sum()
            b[8:] = flipud(b[:7])
        elif method == "KL6":
            s = 10
            a = zeros(s)
            b = zeros(s)
            a[0] =  0.0
            a[1] =  0.39216144400731413927925056
            a[2] =  0.33259913678935943859974864
            a[3] = -0.70624617255763935980996482
            a[4] =  0.08221359629355080023149045
            a[5] =  0.79854399093483008353899777
            a[6:] = flipud(a[1:5])
            b[:s/2] = a[:s/2] + diff(a[:s/2+1])*0.5
            b[s/2:] = flipud(b[:s/2])
        elif method == "KL8":
            s = 18
            a = zeros(s)
            b = zeros(s)
            a[0] =  0.0
            a[1] =  0.13020248308889008087881763
            a[2] =  0.56116298177510838456196441
            a[3] = -0.38947496264484728640807860
            a[4] =  0.15884190655515560089621075
            a[5] = -0.39590389413323757733623154
            a[6] =  0.18453964097831570709183254
            a[7] =  0.25837438768632204729397911
            a[8] =  0.29501172360931029887096624
            a[9] = -0.60550853383003451169892108
            a[10:] = flipud(a[1:9])
            b[:s/2] = a[:s/2] + diff(a[:s/2+1])*0.5
            b[s/2:] = flipud(b[:s/2])
        elif method == "KL10":
            s = 34
            a = zeros(s)
            b = zeros(s)
            a[0]  =  0.0
            a[1]  =  0.09040619368607278492161150
            a[2]  =  0.53591815953030120213784983
            a[3]  =  0.35123257547493978187517736
            a[4]  = -0.31116802097815835426086544
            a[5]  = -0.52556314194263510431065549
            a[6]  =  0.14447909410225247647345695
            a[7]  =  0.02983588609748235818064083
            a[8]  =  0.17786179923739805133592238
            a[9]  =  0.09826906939341637652532377
            a[10] =  0.46179986210411860873242126
            a[11] = -0.33377845599881851314531820
            a[12] =  0.07095684836524793621031152
            a[13] =  0.23666960070126868771909819
            a[14] = -0.49725977950660985445028388
            a[15] = -0.30399616617237257346546356
            a[16] =  0.05246957188100069574521612
            a[17] =  0.44373380805019087955111365
            a[18:] = flipud(a[1:17])
            b[:s/2] = a[:s/2] + diff(a[:s/2+1])*0.5
            b[s/2:] = flipud(b[:s/2])
        else:
            raise NotImplementedError("Unknown method: " + method)

        return a, b


    def intsplit(self, psi1, psi2, a, b, tspan, N, args1=[], args2=[]):
        r"""
        Compute a single, full propagation step by operator splitting.

        :param psi1: First evolution operator :math:`\Psi_a`
        :param psi2: Second evolution operator :math:`\Psi_b`
        :param a: Parameters for evolution with :math:`\Psi_a`
        :param b: Parameters for evolution with :math:`\Psi_b`
        :param tspan: Timespan :math:`t` of a single, full splitting step
        :param N: Number of substeps to perform
        :param args1: Additional optional arguments of :math:`\Psi_a`
        :param args2: Additional optional arguments of :math:`\Psi_b`

        .. note:: The values for ``args1`` and ``args2`` have to be
                  ``list``s even in case of single items.
        """
        s = a.shape[0]
        h = (tspan[1] - tspan[0]) / float(N)

        for k in xrange(N):
            for j in xrange(s):
                psi1(a[j]*h, *args1)
                psi2(b[j]*h, *args2)
