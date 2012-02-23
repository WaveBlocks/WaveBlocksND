"""The WaveBlocks Project

Function for plotting functions of the type f:I^2 -> C
with contour levels having a fixed level distance.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import amax, amin, abs, arange, squeeze
from matplotlib.pyplot import gca


def contourcf(gridx, gridy, values, leveldist=0.02, ax=None):
    """Plot a function :math:`f:I^2 \rightarrow C` with contour levels.
    Put levels in linear and fixed distances on the third axis.

    :param gridx: The grid nodes along the :math:`x` axis of the real domain :math:`R^2`
    :param gridy: The grid nodes along the :math:`y` axis of the real domain :math:`R^2`
    :param values: The values :math:`f(x,y)` of the function.
    :param leveldist: The distance between the contour levels. Default is 0.02.
    """
    # Find the z-values of the contour levels. We keep the
    # z-distance between the levels constant, but not the
    # number of contour levels.
    levels = arange(amin(values), amax(values), leveldist)

    if ax is None:
        ax=gca()

    # Change axes because x.shape = (Nrows, 1) and y.shape = (1, Ncols)
    # is transposed to what we want to plot as x and y axes,
    ax.contour(squeeze(gridy), squeeze(gridx), squeeze(values), levels=levels)
