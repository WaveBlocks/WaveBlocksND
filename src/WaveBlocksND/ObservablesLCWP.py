"""The WaveBlocks Project

Compute some observables like norm, kinetic and potential energy
of linear combinations of general wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import conjugate, transpose, dot, sqrt, array, repeat

from Observables import Observables
from LinearCombinationOfWPs import LinearCombinationOfWPs

__all__ = ["ObservablesLCWP"]


class ObservablesLCWP(Observables):
    r"""This class implements observable computation for linear combinations :math:`\Upsilon`
    of wavepackets :math:`\Psi_j`. There are no assumptions made on the type of the wavepackets
    :math:`\Psi_j` in :math:`\Upsilon := \sum_{j=0}^J c_j \Psi_j`.
    """

    def __init__(self, innerproduct=None):
        r"""Initialize a new :py:class:`ObservablesLCWP` instance for observable computation
        of linear combinations :math:`\Upsilon` of wavepackets :math:`\Psi_j`.

        :param innerproduct: An inner product for computing the integrals. The inner product is
                             used for the computation of brakets :math:`\langle\Psi|\cdot|\Psi\rangle`.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        # A innerproduct to compute the integrals
        if innerproduct is not None:
            self._innerproduct = innerproduct


    def set_innerproduct(self, innerproduct):
        r"""Set the innerproduct.

        :param innerproduct: An inner product for computing the integrals. The inner product is
                             used for the computation of brakets :math:`\langle\Psi|\cdot|\Psi\rangle`.
        :type innerproduct: A :py:class:`InnerProduct` subclass instance.

        .. note:: Make sure to use an inhomogeneous inner product here.
        """
        self._innerproduct = innerproduct


    def overlap_matrix(self, lincomb, component=None):
        r"""Compute the overlap matrix :math:`M_{r,c} = \langle\Upsilon_r|\Upsilon_c\rangle`.

        Note that this function is just a shortcut for calling the inner
        product evaluator directly.

        :param lincomb: The linear combination :math:`\Upsilon`.
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose overlap
                          is calculated. The default value is ``None`` which means to compute
                          overlaps of all :math:`N` components.
        :type component: int or ``None``.
        :return: The matrix :math:`M`.
        """
        OM = self._innerproduct.build_matrix(lincomb, component=component)
        return OM


    def norm(self, lincomb, matrix=None, component=None, summed=False, return_matrix=False):
        r"""Compute the :math:`L^2` norm :math:`\langle\Upsilon|\Upsilon\rangle`
        of a linear combination :math:`\Upsilon` of wavepackets.

        :param lincomb: The linear combination :math:`\Upsilon` of which we compute the norm.
        :type lincomb: A :py:class:`LinearCombinationOfWavepackets` subclass instance.
        :param matrix: The overlap matrix. If ``None`` the matrix is computed internally.
        :type matrix: An ``ndarray`` or ``None`` (default).
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose norm
                          is calculated. The default value is ``None`` which means to compute
                          norms of all :math:`N` components.
        :param return_matrix: Whether to return the overlap matrix used internally.
        :type return_matrix: Boolean, default is ``False``.
        :return: The norm of :math:`\Upsilon` and optionally the overlap matrix :math:`M`.
        """
        if matrix is None:
            OM = self.overlap_matrix(lincomb, component=component)
        else:
            OM = matrix

        if component is not None:
            N = array([1] * lincomb.get_number_packets())
        else:
            N = array([ wp.get_number_components() for wp in lincomb.get_wavepackets() ])

        # Prepare the coefficients in case of multiple components
        c = lincomb.get_coefficients()
        c = repeat(c, N)

        # Compute the norm
        norm = dot(conjugate(transpose(c)), dot(OM, c))
        norm = sqrt(norm)

        # Allow to return the overlap matrix.
        if return_matrix:
            return (norm, OM)
        else:
            return norm


    def kinetic_overlap_matrix(self, lincomb, component=None):
        r"""Compute the kinetic overlap matrix :math:`{M_T}_{r,c} = \langle\Upsilon_r|T|\Upsilon_c\rangle`.

        :param lincomb: The linear combination :math:`\Upsilon`.
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose
                          kinetic energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :return: The matrix :math:`M_T`.
        """
        D = lincomb.get_dimension()
        N = lincomb.get_number_components()

        # TODO: Optimizing this. For large linear combinations, computing
        #       and storing *all* gradient packets is inefficient. Maybe
        #       better go packet by packet.
        grad_lcs = [ LinearCombinationOfWPs(D, N) for d in xrange(D) ]
        for packet in lincomb.get_wavepackets():
            G = packet.get_gradient_operator()
            gradient_wps = G.apply_gradient(packet, as_packet=True)
            # Coefficient c_i=1.0 is wrong but won't be used for building the matrix anyway.
            for d, grad_wp in enumerate(gradient_wps):
                grad_lcs[d].add_wavepacket(grad_wp)

        # Compute the matrices and sum up
        OMT = self._innerproduct.build_matrix(grad_lcs[0], component=component)
        for grad_lc in grad_lcs[1:]:
            OMT += self._innerproduct.build_matrix(grad_lc, component=component)

        return OMT


    def kinetic_energy(self, lincomb, matrix=None, component=None, summed=False, return_matrix=False):
        r"""Compute the kinetic energy :math:`E_{\text{kin}} := \langle\Upsilon|T|\Upsilon\rangle`
        of a linear combination :math:`\Upsilon` of wavepackets.

        :param linbomc: The linear combination :math:`\Upsilon` of which we compute the kinetic energy.
        :type lincomb: A :py:class:`LinearCombinationOfWavepackets` subclass instance.
        :param matrix: The kinetic overlap matrix. If ``None`` the matrix is computed internally.
        :type matrix: An ``ndarray`` or ``None`` (default).
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose
                          kinetic energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param return_matrix: Whether to return the kinetic overlap matrix used internally.
        :type return_matrix: Boolean, default is ``False``.
        :return: The kinetic energy of :math:`\Upsilon` and optionally the kinetic overlap matrix :math:`M_T`.
        """
        if matrix is None:
            OMT = self.kinetic_overlap_matrix(lincomb, component=component)
        else:
            OMT = matrix

        if component is not None:
            N = array([1] * lincomb.get_number_packets())
        else:
            N = array([ wp.get_number_components() for wp in lincomb.get_wavepackets() ])

        # Prepare the coefficients in case of multiple components
        c = lincomb.get_coefficients()
        c = repeat(c, N)

        # Compute the kinetic energy
        ekin = 0.5 * dot(conjugate(transpose(c)), dot(OMT, c))

        # Allow to return the overlap matrix.
        if return_matrix:
            return (ekin, OMT)
        else:
            return ekin


    def potential_overlap_matrix(self, lincomb, potential, component=None):
        r"""Compute the potential overlap matrix :math:`{M_V}_{r,c} = \langle\Upsilon_r|V|\Upsilon_c\rangle`.

        :param lincomb: The linear combination :math:`\Upsilon`.
        :param potential: The potential :math:`V(x)`. (Actually, not the potential object itself
                          but one of its ``V.evaluate_*`` methods.)
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose
                          potential energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :return: The matrix :math:`M_V`.
        """
        OMV = self._innerproduct.build_matrix(lincomb, operator=potential, component=component)
        return OMV


    def potential_energy(self, lincomb, potential, matrix=None, component=None, summed=False, return_matrix=False):
        r"""Compute the potential energy :math:`E_{\text{pot}} := \langle\Upsilon|V|\Upsilon\rangle`.
        of a linear combination :math:`\Upsilon` of wavepackets.

        :param linbomc: The linear combination :math:`\Upsilon` of which we compute the potential energy.
        :type lincomb: A :py:class:`LinearCombinationOfWavepackets` subclass instance.
        :param potential: The potential :math:`V(x)`. (Actually, not the potential object itself
                          but one of its ``V.evaluate_*`` methods.)
        :param matrix: The potential overlap matrix. If ``None`` the matrix is computed internally.
        :type matrix: An ``ndarray`` or ``None`` per default.
        :param component: The index :math:`i` of the components :math:`\Phi_i` whose
                          potential energy we want to compute. If set to ``None`` the
                          computation is performed for all :math:`N` components.
        :type component: Integer or ``None``.
        :param return_matrix: Whether to return the potential overlap matrix used internally.
        :type return_matrix: Boolean, default is ``False``.
        :return: The potential energy of :math:`\Upsilon` and optionally the potential overlap matrix :math:`M_V`.
        """
        if matrix is None:
            OMV = self.potential_overlap_matrix(lincomb, potential, component=component)
        else:
            OMV = matrix

        if component is not None:
            N = array([1] * lincomb.get_number_packets())
        else:
            N = array([ wp.get_number_components() for wp in lincomb.get_wavepackets() ])

        # Prepare the coefficients in case of multiple components
        c = lincomb.get_coefficients()
        c = repeat(c, N)

        # Compute the potential energy
        epot = dot(conjugate(transpose(c)), dot(OMV, c))

        # Allow to return the overlap matrix.
        if return_matrix:
            return (epot, OMV)
        else:
            return epot
