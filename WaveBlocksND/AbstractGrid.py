"""The WaveBlocks Project

This file contains the abstract class for representing grids.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

__all__ = ["AbstractGrid"]


class AbstractGrid(object):
    """This class is an abstract interface to grids in general.
    """

    def __init__(self):
        """
        :raise: :py:class:`NotImplementedError` Abstract interface.
        """
        raise NotImplementedError("'Grid' is an abstract interface.")


    def get_dimension(self):
        """Return the dimension :math:`D` of the grid.
        """
        return self._dimension
