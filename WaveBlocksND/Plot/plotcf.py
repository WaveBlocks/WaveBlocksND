"""The WaveBlocks Project

Plot function for plotting functions of the type f:R -> C
with abs(f) as y-value and phase(f) as color code.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import array
from matplotlib.collections import LineCollection
from matplotlib.pyplot import gca

from WaveBlocksND.Plot.color_map import color_map


def plotcf(grid, phase, modulus, darken=None, axes=None, linestylep="solid", linewidthp=1, color="k", **kwargs):
    r"""Plot the modulus of a complex valued function
    :math:`f:\mathbb{R} \rightarrow \mathbb{C}`
    together with its phase in a color coded fashion.

    :param grid: The grid nodes of the real domain R
    :param phase: The phase of the complex domain result f(grid)
    :param modulus: The modulus of the complex domain result f(grid)
    :param darken: Whether to take into account the modulus of the data to darken colors.
    :param axes: The axes instance used for plotting.
    :param linestylep: The line style of the phase curve.
    :param linewidthp: The line width of the phase curve.
    :param color: The color of the phase curve.
    """
    # Color mapping
    rgb_colors = color_map(modulus, phase=phase, modulus=modulus, darken=darken)

    # Put all the vertical line into a collection
    segments = [array([[node, 0], [node, value]]) for node, value in zip(grid, modulus)]
    line_segments = LineCollection(segments)

    # Set some properties of the lines
    rgb_colors = line_segments.to_rgba(rgb_colors)
    line_segments.set_color(rgb_colors[0])
    line_segments.set_linestyle(linestylep)
    line_segments.set_linewidth(linewidthp)

    # Plot to the given axis instance or retrieve the current one
    if axes is None:
        axes = gca()

    # Plot the phase
    axes.add_collection(line_segments)
    # Plot the modulus
    axes.plot(grid, modulus, color=color, **kwargs)
