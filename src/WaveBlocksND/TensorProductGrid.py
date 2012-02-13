"""The WaveBlocks Project

This file contains a class for representing
dense regular tensor product grids.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import operator
from numpy import mgrid, ogrid, atleast_1d

from DenseGrid import DenseGrid

__all__ = ["TensorProductGrid"]


class TensorProductGrid(DenseGrid):
    """This class represents a dense tensor product grid. It can
    have an arbitrary dimension :math:`D`. The grid nodes are
    enclosed in a hypercubic bounding box. This box can have
    different limits :math:`min_i`, :math:`max_i` along each
    axis :math:`x_i`. In each of these intervalls we place
    :math:`N_i` grid nodes. Note that the point :math:`max_i`
    is not part of the grid. The grid interior is build as the
    tensor product of all the grid nodes along all the axes.
    """

    def __init__(self, parameters):
        """Construct a tnesor product grid instace.

        :param parameters: The simulation parameters. There must be
                           values for the keys `dimension`, `grid_limits`
                           and `grid_number_nodes`.
        :return: A :py:class:`TensorProductGrid` instance.
        """

        # TODO: Improve for one|multiple values
        #       This comes with the Blocks factory?

        # Regular grid spacing
        self.is_regular = True

        # The dimension of the grid
        self.dimension = parameters["dimension"]

        # The limits of the bounding box of the grid
        self.limits = parameters["grid_limits"]
        # format: [(min_0,max_0), ..., (min_D,max_D)]

        # The number of grid nodes along each axis
        self.number_nodes = parameters["grid_number_nodes"]
        # format: [N_1, ..., N_D]

        # Compute the grid spacings along each axis
        self.meshwidth = []
        # format: [h_1, ...,, h_D]
        for lims, Ni in zip(self.limits, self.number_nodes):
            self.meshwidth.append(abs(lims[1]-lims[0])/(1.0*Ni))

        # Cached values
        self._gridaxes = None
        self._gridnodes = None


    def get_limits(self, axes=None):
        """Returns the limits of the bounding box.

        :param axes: The axes for which we want to get the limits.
        :type axes: A single integer or a list of integers. If set
                    to `None` (default) we return the data for all axes.
        :return: A list of :math:`(min_i, max_i)` tuples or a single tuple.
        """
        if axes is None:
            axes = range(self.dimension)

        limits = [ self.limits[i] for i in atleast_1d(axes) ]

        if len(limits) == 1:
            return limits[0]
        return limits


    def get_meshwidth(self, axes=None):
        """Returns the meshwidth of the grid.

        :param axes: The axes for which we want to get the meshwidth.
        :type axes: A single integer or a list of integers. If set
                    to `None` (default) we return the data for all axes.
        :return: A list of :math:`h_i` values or a single value.
        """
        if axes is None:
            axes = xrange(self.dimension)

        meshwidth = [ self.meshwidth[i] for i in atleast_1d(axes) ]

        if len(meshwidth) == 1:
            return meshwidth[0]
        return meshwidth


    def get_number_nodes(self, axes=None, overall=True):
        """Returns the number of grid nodes along a set of axes.

        :param axes: The axes for which we want to get the number of nodes.
        :type axes: A single integer or a list of integers. If set
                    to `None` (default) we return the data for all axes.
        :param overall: Compute the product :math:`\prod_i N_i` of the
                        number :math:`N_i` of grid nodes along each axis
                        :math:`i` specified.
        :type overall: Boolean, default is `True`
        :return: A list of :math:`h_i` values or a single value.
        """
        if axes is None:
            axes = xrange(self.dimension)

        values = [ self.number_nodes[i] for i in atleast_1d(axes) ]

        if overall is True:
            return reduce(operator.mul, values)
        elif len(values) == 1:
            return values[0]
        else:
            return values


    def _build_slicers(self):
        # Helper routine to build the necessary slicing
        # objects used for constructing the grid nodes.
        slicers = [ slice(lims[0], lims[1], step) for lims, step in zip(self.limits, self.meshwidth) ]
        return slicers


    def _compute_grid_axes(self):
        # Helper routine which computes the one-dimensional
        # grids along all axes and caches the result. Each
        # grid is a D-dimensional ndarray of correct shape.
        S = self._build_slicers()
        self._gridaxes = ogrid[S]


    def _compute_grid_full(self):
        # Helper routine which computes the full set of
        # tensor product grid nodes and caches the result.
        # The result is a (D, product(N_i)) shaped ndarray.
        S = self._build_slicers()
        grid =  mgrid[S]
        # TODO: Consider storing a (D, N_1, ..., N_D) shaped version
        self._gridnodes = grid.reshape((self.dimension, self.get_number_nodes()))


    def get_axes(self, axes=None):
        """Returns the one-dimensional grids along the axes.

        :param axes: The axes for which we want to get the grid.
        :type axes: A single integer or a list of integers. If set
                    to `None` (default) we return the data for all axes.
        :return: A list of ndarrays, each having a shape of :math:`(1,...,N_i,...,1)`.
                 We return a list even if it contains just a single element.
        """
        if self._gridaxes is None:
            self._compute_grid_axes()

        if axes is None:
            axes = xrange(self.dimension)
        axes = atleast_1d(axes)

        return [ self._gridaxes[i] for i in axes ]


    def get_nodes(self):
        """Returns all grid nodes of the full tensor product grid.

        :return: An ndarrays of shape :math:`(D, \prod_i^D N_i)`.
        """
        if self._gridnodes is None:
            self._compute_grid_full()

        return self._gridnodes
