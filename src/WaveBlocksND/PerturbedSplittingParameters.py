"""The WaveBlocks Project

This file contains data to build several closely
related splitting methods for perturbed operators.

@author: V. Gradinaru, R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, floating

__all__ = ["PerturbedSplittingParameters"]


class PerturbedSplittingParameters(object):


    def build(self, method):
        r"""
        :param method: A string specifying the method for time integration.
        :return: Two arrays :math:`a` and :math:`b`.

        ====== ====== ========= =========
        Method Order  Authors   Reference
        ====== ====== ========= =========
        L42    (4,2)  McLachlan [1]_ page 6
        L62    (6,2)  McLachlan [1]_ page 6
        L82    (8,2)  McLachlan [1]_ page 6
        L10    (10,2) McLachlan [1]_ page 6
        L84    (8,4)  McLachlan [1]_ page 8
        ====== ====== ========= =========

        .. [1] R.I. McLachlan, "Composition methods in the presence of small parameters",
               BIT Numerical Mathematics, Volume 35, Issue 2, (1995) 258-268.
               There is also a 2003 version which was used for the implementation.
        """
        if method == "L42":
            # Pattern ABA and m = s = 2
            m = 2
            a = zeros(m+1)
            b = zeros(m)
            #
            a[0] = 0.21132486540518713
            a[1] = 0.57735026918962584
            a[2] = 0.21132486540518713
            #
            b[0] = 0.5
            b[1] = 0.5
        elif method == "L62":
            # Pattern ABA and m = s = 3
            m = 3
            a = zeros(m+1)
            b = zeros(m)
            #
            a[0] = 0.1127016653792583
            a[1] = 0.3872983346207417
            a[2] = 0.3872983346207417
            a[3] = 0.1127016653792583
            #
            b[0] = 0.2777777777777778
            b[1] = 0.4444444444444444
            b[2] = 0.2777777777777778
        elif method == "L82":
            # Pattern ABA and m = s = 4
            m = 4
            a = zeros(m+1)
            b = zeros(m)
            #
            a[0] = 0.069431844202973714
            a[1] = 0.26057763400459816
            a[2] = 0.33998104358485626
            a[3] = 0.26057763400459816
            a[4] = 0.069431844202973714
            #
            b[0] = 0.17392742256872692
            b[1] = 0.3260725774312731
            b[2] = 0.3260725774312731
            b[3] = 0.17392742256872692
        elif method == "L102":
            # Pattern ABA and m = s = 5
            m = 5
            a = zeros(m+1)
            b = zeros(m)
            #
            a[0] = 0.046910077030668018
            a[1] = 0.18385526791649043
            a[2] = 0.26923465505284155
            a[3] = 0.26923465505284155
            a[4] = 0.18385526791649043
            a[5] = 0.046910077030668018
            #
            b[0] = 0.11846344252809454
            b[1] = 0.23931433524968324
            b[2] = 0.28444444444444444
            b[3] = 0.23931433524968324
            b[4] = 0.11846344252809454
        elif method == "L84":
            # Pattern ABA and m = 5
            m = 5
            a = zeros(m+1)
            b = zeros(m)
            #
            a[0] =  0.07534696026989288842
            a[1] =  0.51791685468825678230
            a[2] = -0.09326381495814967072
            a[3] = -0.09326381495814967072
            a[4] =  0.51791685468825678230
            a[5] =  0.07534696026989288842
            #
            b[0] =  0.19022593937367661925
            b[1] =  0.84652407044352625706
            b[2] = -1.07350001963440575260
            b[3] =  0.84652407044352625706
            b[4] =  0.19022593937367661925
        else:
            raise NotImplementedError("Unknown method: " + method)

        return a, b


    def intsplit(self, psia, psib, a, b, tspan, N, argsa=[], argsb=[]):
        r"""
        Compute a single, full propagation step by operator splitting.

        :param psia: First evolution operator :math:`\Psi_a`
        :param psib: Second evolution operator :math:`\Psi_b`
        :param a: Parameters for evolution with :math:`\Psi_a`
        :param b: Parameters for evolution with :math:`\Psi_b`
        :param tspan: Timespan :math:`t` of a single, full splitting step
        :param N: Number of substeps to perform
        :param argsa: Additional optional arguments of :math:`\Psi_a`
        :param argsb: Additional optional arguments of :math:`\Psi_b`

        .. note:: The values for ``argsa`` and ``argsb`` have to be
                  ``list``s even in case of single items.
        """
        psi = (psia, psib)
        args = (argsa, argsb)
        h = (tspan[1] - tspan[0]) / float(N)
        sa = a.shape[0]
        sb = b.shape[0]
        c = zeros(sa+sb, dtype=floating)
        c[::2] = a
        c[1::2] = b

        for k in xrange(N):
            for j in xrange(sa+sb):
                psi[j%2](c[j]*h, *args[j%2])
