"""The WaveBlocks Project

Function for stem-plotting functions of the type f:I -> C
with abs(f) as y-value and phase(f) as color code.
This function makes a stem plot.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import array, zeros, real
from matplotlib.collections import LineCollection
from matplotlib.pyplot import gca

from color_map import color_map


def stemcf(grid, phase, modulus, darken=False, axes=None, linestylep="solid", linewidthp=2, color=None, markerp="o", **kwargs):
    r"""Stemplot the modulus of a complex valued function :math:`f:I -> \mathbb{C}` together with its phase in a color coded fashion.
    Additional keyword arguments are passed to the plot function.

    :param grid: The grid nodes of the real domain grid :math:`\Gamma`
    :param phase: The phase of the complex domain result :math:`f(\Gamma)`
    :param modulus: The modulus of the complex domain result :math:`f(\Gamma)`
    :param darken: Whether to take into account the modulus of the data to darken colors.
    :param axes: The axes instance used for plotting.
    :param linestylep: The line style of the phase curve.
    :param linewidthp: The line width of the phase curve.
    :param color: The color of the stemmed markers.
    :param markerp: The shape of the stemmed markers.
    """
    # Color mapping
    rgb_colors = color_map(grid, phase=phase, modulus=modulus, darken=darken)

    # Put all the vertical line into a collection
    segments = [ array([[node,0], [node,value]]) for node, value in zip(grid, modulus) ]
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
    if color is None:
        # Scatter has a problem with complex data type, make sure values are purely real
        axes.scatter(grid, real(modulus), c=rgb_colors[0], **kwargs)
    else:
        axes.plot(grid, modulus, linestyle="", marker=markerp, color=color, **kwargs)
    # Plot the ground line
    axes.plot(grid, zeros(grid.shape), linestyle=linestylep, color="k", **kwargs)
