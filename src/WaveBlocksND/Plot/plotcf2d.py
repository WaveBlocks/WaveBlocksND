"""The WaveBlocks Project

Function for plotting complex valued functions
of two real variables with the values encoded
by the usual color code.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

from numpy import meshgrid, real, where
from matplotlib.pyplot import gca

from .color_map import color_map


def plotcf2d(x, y, z, darken=None, axes=None, limits=None, **kwargs):
    r"""Plot complex valued functions :math:`\mathbb{R}^2 \rightarrow \mathbb{C}`
    with the usual color code.

    :param x: The :math:`x` values.
    :param x: The :math:`y` values.
    :param z: The values :math:`z = f(x,y)`.
    :param darken: How strong to take into account the modulus of the data to darken colors.
                   Values with :math:`|z| = R` will get fully saturated colors
                   while :math:`|z| = 0` is black and :math:`|z| \rightarrow \infty`
                   get whiter and whiter.
    :type darken: Float or ``None`` to disable darkening of colors. Default is :math:`R = 1.0`.
    :param axes: The axes instance used for plotting.
    """
    if limits is None:
        xmin = real(x).min()
        xmax = real(x).max()
        ymin = real(y).min()
        ymax = real(y).max()
        extent = [xmin, xmax, ymin, ymax]
    else:
        xmin = limits[0]
        xmax = limits[1]
        ymin = limits[2]
        ymax = limits[3]
        extent = [xmin, xmax, ymin, ymax]

    kw = {'extent': extent,
          'origin': 'lower',
          'interpolation': 'nearest',
          'aspect': 'equal'}
    kw.update(kwargs)

    # Plot to the given axis instance or retrieve the current one
    if axes is None:
        axes = gca()

    # Region to cut out
    x = x.reshape(1,-1)
    y = y.reshape(-1,1)
    i = where((xmin <= x) & (x <= xmax))[1]
    j = where((ymin <= y) & (y <= ymax))[0]
    I, J = meshgrid(i,j)

    # Color code and plot the data
    axes.imshow(color_map(z[I,J], darken=darken), **kw)
