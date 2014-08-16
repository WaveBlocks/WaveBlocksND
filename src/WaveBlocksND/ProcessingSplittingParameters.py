"""The WaveBlocks Project

This file contains data to build several closely
related processing splitting methods.

@author: V. Gradinaru, R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros

__all__ = ["ProcessingSplittingParameters"]


class ProcessingSplittingParameters(object):


    def build(self, method):
        r"""
        :param method: A string specifying the method for time integration.
        :return: Two arrays :math:`a` and :math:`b`.

        ====== ======= ================= =========
        Method Order   Authors           Reference
        ====== ======= ================= =========
        BCR764 (7,6,4) Blanes/Casas/Ros  [1]_ table (iv)
        ====== ======= ================= =========

        .. [1] Blanes, Casas, Ros 2000, table IV
        """
        if method == "BCR764":
            # Kernel Pattern ABA
            # Exchanged a and b compared to the paper in
            # order to have consistency with other code
            a = zeros(4)
            a[1] = 1.5171479707207228
            a[2] = 1.0 - 2.0 * a[1]
            a[3] = 1.0 * a[1]
            b = zeros(4)
            b[0] = 0.5600879810924619
            b[1] = 0.5 - b[0]
            b[2] = 1.0 * b[1]
            b[3] = 1.0 * b[0]
            # Pre/Post-Processor
            z = zeros(6)
            z[0] = -0.3346222298730
            z[1] = 1.097567990732164
            z[2] = -1.038088746096783
            z[3] = 0.6234776317921379
            z[4] = -1.102753206303191
            z[5] = -0.0141183222088869
            y = zeros(6)
            y[0] = -1.621810118086801
            y[1] = 0.0061709468110142
            y[2] = 0.8348493592472594
            y[3] = -0.0511253369989315
            y[4] = 0.5633782670698199
            y[5] = -0.5
        else:
            raise NotImplementedError("Unknown method: " + method)

        return a, b, y, z


    def intprepsplit(self, psi1, psi2, a, b, y, z, tspan, N, args1=[], args2=[]):
        r"""
        Compute a single, full propagation step by processing operator splitting.

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
        p = y.shape[0]
        h = (tspan[1] - tspan[0]) / float(N)

        # Preprocessor
        for j in xrange(p):
            psi1(-z[j]*h, *args1)
            psi2(-y[j]*h, *args2)

        # Kernel
        for k in xrange(N):
            for j in xrange(s):
                psi1(a[j]*h, *args1)
                psi2(b[j]*h, *args2)

        # Postprocessor
        for j in xrange(p-1, -1, -1):
            psi1(y[j]*h, *args1)
            psi2(z[j]*h, *args2)
