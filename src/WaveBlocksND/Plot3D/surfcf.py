r"""The WaveBlocks Project

Function for plotting functions of the type :math:`f:I^2 -> C`
with :math:`|f|` as y-value and :math:`\arg(f)` as color code.
This function makes a 3D surface plot.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import linspace, pi, squeeze, ones, real, fmod
from matplotlib.colors import hsv_to_rgb
from mayavi import mlab


def compute_color_map():
    """Compute a default QM colormap which can be used as mayavi/vtk lookup table.
    """
    k = linspace(-pi, pi, 256, endpoint=True)
    hsv_colors = ones((1, k.shape[0], 3))
    hsv_colors[:,:, 0] = 0.5*fmod(k+2*pi,2*pi)/pi
    return 255*squeeze(hsv_to_rgb(hsv_colors))


def surfcf(gridx, gridy, phase, modulus, colormap=None):
    r"""Plot the modulus of a complex valued function :math:`f:R^2 -> C`
    together with its phase in a color coded fashion.

    :param gridx: The grid nodes along the :math:`x` axis of the real domain :math:`R^2`
    :param gridy: The grid nodes along the :math:`y` axis of the real domain :math:`R^2`
    :param phase: The phase of the complex domain result f(grid)
    :param modulus: The modulus of the complex domain result f(grid)
    :param colormap: The colormap to use, if none is given, compute the 'default' QM colormap.
    """
    if colormap is None:
        colormap = compute_color_map()

    # The real(.) is necessary just to get an array with dtype real
    mesh = mlab.mesh(gridx, gridy, real(modulus), scalars=phase)

    # Set the custom color map
    mesh.module_manager.scalar_lut_manager.use_default_range = False
    mesh.module_manager.scalar_lut_manager.data_range = [-pi, pi]
    lut = mesh.module_manager.scalar_lut_manager.lut.table.to_array()
    lut[:,0:3] = colormap.copy()
    mesh.module_manager.scalar_lut_manager.lut.table = lut

    # Update the figure
    mlab.draw()

    return mesh
