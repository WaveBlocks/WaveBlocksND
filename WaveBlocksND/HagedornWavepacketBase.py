"""The WaveBlocks Project

This file contains the basic interface for general wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import array, complexfloating, cumsum, vsplit, vstack, zeros, conjugate
from scipy import sqrt
from scipy.linalg import norm

from WaveBlocksND.Wavepacket import Wavepacket
from WaveBlocksND.GradientHAWP import GradientHAWP
from functools import reduce

__all__ = ["HagedornWavepacketBase"]


class HagedornWavepacketBase(Wavepacket):
    r"""This class implements the abstract :py:class:`Wavepacket` interface
    and contains code common to all types of Hagedorn wavepackets.
    """

    def __init__(self, parameters):
        r"""Initialize a wavepacket object that represents :math:`\Psi`.
        """
        raise NotImplementedError("'HagedornWavepacketBase' should not be instantiated.")


    # We can handle basis shapes here as the logic is the same for
    # homogeneous and inhomogeneous Hagedorn wavepackets.

    def _resize_coefficient_storage(self, component, bs_old, bs_new):
        r"""
        """
        bso = bs_old.get_basis_size()
        bsn = bs_new.get_basis_size()

        # Find the intersection of K and K'
        # Optimization: iterate over smaller set
        if bso <= bsn:
            insec = [k for k in bs_old if k in bs_new]
        elif bso > bsn:
            insec = [k for k in bs_new if k in bs_old]
        # TODO: Consider making this part of the BasisShape interface
        # TODO: Consider implementing set operations for basis shapes

        # Construct the index mapping
        i = array([bs_old[k] for k in insec])
        j = array([bs_new[k] for k in insec])

        # Copy over the data
        cnew = zeros((bsn, 1), dtype=complexfloating)
        cnew[j] = self._coefficients[component][i]
        self._coefficients[component] = cnew


    def get_basis_shapes(self, *, component=None):
        r"""Retrieve the basis shapes :math:`\mathfrak{K}_i` for each component :math:`i`.

        :param component: The component :math:`i` whose basis shape we request. (Default is
                          ``None`` which means to return the basis shapes for all components.
        :type component: int
        :return: The basis shape for an individual component or a list with all shapes.
        """
        if component is not None:
            return self._basis_shapes[component]
        else:
            return tuple(self._basis_shapes)


    def set_basis_shapes(self, basis_shape, *, component=None):
        r"""Set the basis shape :math:`\mathfrak{K}` of a given component or for all components.

        :param basis_shape: The basis shape for an individual component or a list with all :math:`N` shapes.
        :type basis_shape: A subclass of :py:class:`BasisShape`.
        :param component: The component :math:`i` whose basis shape we want to set. (Default is
                          ``None`` which means to set the basis shapes for all components.
        :type component: int
        """
        if component is not None:
            # Check for valid input basis shape
            if component not in range(self._number_components):
                raise ValueError("Invalid component index " + str(component))

            # Adapt the coefficient storage vectors
            self._resize_coefficient_storage(component, self._basis_shapes[component], basis_shape)
            # Set the new basis shape for the given component
            self._basis_shapes[component] = basis_shape
        else:
            # Check for valid input basis shape
            if not len(basis_shape) == self._number_components:
                raise ValueError("Number of basis shape(s) does not match to number of components.")

            for index, bsnew in enumerate(basis_shape):
                # Adapt the coefficient storage vectors
                self._resize_coefficient_storage(index, self._basis_shapes[index], bsnew)
                # Set the new basis shape for the given component
                self._basis_shapes[index] = bsnew

        # And update the caches information
        self._basis_sizes = [bs.get_basis_size() for bs in self._basis_shapes]


    # We can handle coefficient set manipulation here as the logic is
    # the same for homogeneous and inhomogeneous Hagedorn wavepackets.


    def set_coefficient(self, component, index, value):
        r"""Set a single coefficient :math:`c^i_k` of the specified component :math:`\Phi_i`
        of :math:`\Psi`.

        :param component: The index :math:`i` of the component :math:`\Phi_i` we want to update.
        :type components: int
        :param index: The multi-index :math:`k` of the coefficient :math:`c^i_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :param value: The new value of the coefficient :math:`c^i_k`.
        :raise: :py:class:`ValueError` For invalid indices :math:`i` or :math:`k`.
        """
        if component > self._number_components - 1 or component < 0:
            raise ValueError("There is no component with index {}.".format(component))
        if index not in self._basis_shapes[component]:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shapes[component][index]
        self._coefficients[component][key] = value


    def get_coefficient(self, component, index):
        r"""Retrieve a single coefficient :math:`c^i_k` of the specified component :math:`\Phi_i`
        of :math:`\Psi`.

        :param component: The index :math:`i` of the component :math:`\Phi_i` we want to update.
        :type components: int
        :param index: The multi-index :math:`k` of the coefficient :math:`c^i_k` we want to update.
        :type index: A tuple of :math:`D` integers.
        :return: A single complex number.
        :raise: :py:class:`ValueError` For invalid indices :math:`i` or :math:`k`.
        """
        if component > self._number_components - 1 or component < 0:
            raise ValueError("There is no component with index {}.".format(component))
        if index not in self._basis_shapes[component]:
            raise ValueError("There is no basis function with multi-index {}.".format(index))

        # Apply linear order mapping here
        key = self._basis_shapes[component][index]
        return self._coefficients[component][key]


    def set_coefficients(self, values, *, component=None):
        r"""Update all the coefficients :math:`c` of :math:`\Psi` or update
        the coefficients :math:`c^i` of the components :math:`\Phi_i` only.

        Note: this method copies the data arrays.

        :param values: The new values of the coefficients :math:`c^i` of :math:`\Phi_i`.
        :type values: An ndarray of suitable shape or a list of ndarrays.
        :param component: The index :math:`i` of the component we want to update with new coefficients.
        :type component: int (Default is ``None`` meaning all)
        :raise: :py:class:`ValueError` For invalid component indices :math:`i`.
        """
        if component is None:
            if len(values) != self._number_components:
                raise ValueError("Too less or too many data provided.")

            for index, value in enumerate(values):
                bs = self._basis_sizes[index]
                self._coefficients[index] = value.copy().reshape((bs, 1))
        else:
            if component > self._number_components - 1 or component < 0:
                raise ValueError("There is no component with index {}.".format(component))

            bs = self._basis_sizes[component]
            self._coefficients[component] = values.copy().reshape((bs, 1))


    def get_coefficients(self, *, component=None):
        r"""Returns the coefficients :math:`c^i` for some component :math:`\Phi_i` of
        :math:`\Psi` or all the coefficients :math:`c` of all components.

        Note: this method copies the data arrays.

        :param component: The index :math:`i` of the component we want to retrieve.
        :type component: int (Default is ``None`` meaning all)
        :return: A single ndarray with the coefficients of the given component or
                 a list containing the ndarrays for each component. Each ndarray
                 is two-dimensional with a shape of :math:`(|\mathfrak{K}_i|, 1)`.
        :raise: :py:class:`ValueError` For invalid component indices :math:`i`.
        """
        if component is None:
            return [item.copy() for item in self._coefficients]
        else:
            if component > self._number_components - 1 or component < 0:
                raise ValueError("There is no component with index {}.".format(component))

            return self._coefficients[component].copy()


    def get_coefficient_vector(self, *, component=None):
        r"""Retrieve the coefficients for all components :math:`\Phi_i` simultaneously.

        .. warning:: This function does *not* copy the input data!
                     This is for efficiency as this routine is used in the innermost loops.

        :param component: The component :math:`i` whose coefficients we request. (Default is
                          ``None`` which means to return the coefficients for all components.
        :type component: int
        :return: The coefficients :math:`c^i` of all components
                 :math:`\Phi_i` stacked into a single long column vector.
        """
        if component is None:
            return vstack(self._coefficients)
        else:
            return self._coefficients[component]


    def set_coefficient_vector(self, vector):
        """Set the coefficients for all components :math:`\Phi_i` simultaneously.

        .. warning:: This function does *not* copy the input data!
                     This is for efficiency as this routine is used in the innermost loops.

        :param vector: The coefficients of all components as a single long column vector.
        :type vector: A two-dimensional ndarray of appropriate shape.
        """
        # Compute the partition of the block-vector from the basis sizes
        partition = cumsum(self._basis_sizes)[:-1]

        # Split the block-vector with the given partition and assign
        self._coefficients = vsplit(vector, partition)


    def get_eps(self):
        r"""Retrieve the semi-classical scaling parameter :math:`\varepsilon` of the wavepacket.

        :return: The value of :math:`\varepsilon`.
        """
        return self._eps


    # We can compute the norms the same way for homogeneous and inhomogeneous Hagedorn wavepackets.


    def norm(self, component=None, summed=False):
        r"""Calculate the :math:`L^2` norm :math:`\langle\Psi|\Psi\rangle` of the wavepacket :math:`\Psi`.

        :param component: The index :math:`i` of the component :math:`\Phi_i` whose norm is calculated.
                          The default value is ``None`` which means to compute the norms of all :math:`N` components.
        :type component: int or ``None``.
        :param summed: Whether to sum up the norms :math:`\langle\Phi_i|\Phi_i\rangle` of the
                       individual components :math:`\Phi_i`.
        :type summed: Boolean, default is ``False``.
        :type summed: Boolean, default is ``False``.
        :return: The norm of :math:`\Psi` or the norm of :math:`\Phi_i` or a list with the :math:`N`
                 norms of all components. Depending on the values of ``component`` and ``summed``.
        """
        if component is not None:
            result = norm(self._coefficients[component])
        else:
            result = [norm(item) for item in self._coefficients]

            if summed is True:
                result = reduce(lambda x, y: x + conjugate(y) * y, result, 0.0)
                result = sqrt(result)

        return result


    # A wavepacket knows how to compute gradients
    # TODO: Consider moving this inside the general codata framework


    def get_gradient_operator(self):
        r"""Return the :py:class:`Gradient` subclass suitable for
        computing gradients of this wavepacket.

        :return: A :py:class:`GradientHAWP` instance.
        """
        return GradientHAWP()


    # A wavepacket knows how to compute inner products
    # TODO: Consider moving this inside the general codata framework
    # TODO: Rethink if wavepackets should contain a QR


    def set_innerproduct(self, innerproduct):
        """Set the :py:class:`InnerProduct` subclass instance used for computing
        inner products and evaluating brakets.

        :param innerproduct: The new :py:class:`InnerProduct` subclass instance.
        """
        self._IP = innerproduct


    def get_innerproduct(self):
        """Return the :py:class:`InnerProduct` subclass instance used computing
        inner products and evaluating brakets.

        :return: The current :py:class:`InnerProduct` subclass instance.
        """
        return self._IP
