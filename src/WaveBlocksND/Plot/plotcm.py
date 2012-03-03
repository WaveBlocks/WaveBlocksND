"""The WaveBlocks Project

Function for plotting complex matrices with the
phase of the entries encoded into the usual color code.

@author: R. Bourquin
@copyright: Copyright (C) 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

from matplotlib.pyplot import gca
from matplotlib import ticker

from color_map import color_map


def plotcm(matrix, phase=None, modulus=None, darken=True, axes=None, **kwargs):
    """Plot complex matrices with the phase of the entries encoded into the usual color code.
    @param matrix: The matrix data.
    @param phase: The phase of the entries, if not given they are computed.
    @param modulus: The modulus of the entries, if not given they are computed.
    @keyword darken: Whether to take into account the modulus of the data to darken colors.
    @keyword axes: The axes instance used for plotting.
    @note: Additional keyword arguments are passe to the plot function.
    """
    # TODO: Allow to specify axes

    nr, nc = matrix.shape
    extent = [-0.5, nc-0.5, nr-0.5, -0.5]

    kw = {'extent': extent,
          'origin': 'upper',
          'interpolation': 'nearest',
          'aspect': 'equal'}
    kw.update(kwargs)

    # Plot to the given axis instance or retrieve the current one
    if axes is None:
        axes = gca()

    # Color code and plot the data matrix
    axes.imshow(color_map(matrix, phase=phase, modulus=modulus, darken=darken), **kw)

    axes.xaxis.tick_top()
    axes.xaxis.set_ticks_position('both')
    axes.xaxis.set_major_locator(ticker.MaxNLocator(nbins=9, steps=[1, 2, 5, 10], integer=True))
    axes.yaxis.set_major_locator(ticker.MaxNLocator(nbins=9, steps=[1, 2, 5, 10], integer=True))
