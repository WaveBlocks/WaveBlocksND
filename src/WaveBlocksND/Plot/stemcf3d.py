"""The WaveBlocks Project

Function for stem-plotting functions of the type f:IxI -> C
with abs(f) as z-value and phase(f) as color code.
This function makes a three dimensional stem plot.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import real, squeeze
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.pyplot import gcf

from color_map import color_map


def stemcf3d(gridu, gridv, phase, modulus, darken=None, fig=None, markerp="o", **kwargs):
    r"""Stemplot the modulus of a complex valued function :math:`f:I\times I -> \mathbb{C}` together with its
    phase in a color coded fashion. Additional keyword arguments are passed to the plot function.

    :param gridu: The x components of the grid nodes of the real domain grid :math:`\Gamma`
    :param gridv: The y components of the grid nodes of the real domain grid :math:`\Gamma`
    :param phase: The phase of the complex domain result :math:`f(\Gamma)`
    :param modulus: The modulus of the complex domain result :math:`f(\Gamma)`
    :param darken: Whether to take into account the modulus of the data to darken colors.
    :param fig: The figure instance used for plotting.
    :param markerp: The shape of the stemmed markers.
    """
    # Color mapping
    rgb_colors = squeeze(color_map(gridv, phase=phase, modulus=modulus, darken=darken))

    # Plot to the given axis instance or retrieve the current one
    if fig is None:
        fig = gcf()

    axes = fig.add_subplot(1, 1, 1, projection='3d')

    for ui, vi, wi, col in zip(gridu, gridv, modulus, rgb_colors):
        line = art3d.Line3D(*zip((ui, vi, 0), (ui, vi, wi)), marker=markerp, markevery=(1, 1), color=col)
        axes.add_line(line)

    axes.set_xlim3d(real(gridu).min(), real(gridu).max())
    axes.set_ylim3d(real(gridv).min(), real(gridv).max())
    axes.set_zlim3d(real(modulus).min(), real(modulus).max())
