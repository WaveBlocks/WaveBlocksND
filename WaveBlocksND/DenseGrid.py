"""The WaveBlocks Project

This file contains the abstract class for representing dense grids.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND.AbstractGrid import AbstractGrid

__all__ = ["DenseGrid"]


class DenseGrid(AbstractGrid):
    """This class is an abstract interface to dense grids in general.
    """

    def __init__(self):
        """
        :raise: :py:class:`NotImplementedError` Abstract interface.
        """
        raise NotImplementedError("'DenseGrid' is an abstract interface.")


    def is_regular(self):
        """Answers the question if the grid spacing is regular.
        (It can still be different along each axis!)
        """
        return self._is_regular
