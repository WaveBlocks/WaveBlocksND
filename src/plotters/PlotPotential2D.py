"""The WaveBlocks Project

Plot the eigenvalues (energy levels) of the potential.
This script is only for two-dimensional potentials.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from numpy import real
from mayavi import mlab

from WaveBlocksND import PotentialFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager


def plot_potential(grid, potential, along_axes=False, interactive=False, size=(800,700)):
    # The Grid
    u, v = grid.get_nodes(split=True, flat=False)
    u = real(u)
    v = real(v)

    # TODO: Consider taking every n-th gridpoint only

    # Create potential and evaluate eigenvalues
    potew = potential.evaluate_eigenvalues_at(grid)

    # Plot the energy surfaces of the potential
    fig = mlab.figure(size=size)

    for level in potew:
        mlab.surf(u, v, level)

    fig.scene.parallel_projection = True
    fig.scene.isometric_view()
    fig.scene.show_axes = True

    mlab.savefig("potential_3D_view.png")

    # Parallele views
    if along_axes is True:
        fig.scene.x_minus_view()
        mlab.savefig("potential_xm_view.png")

        fig.scene.x_plus_view()
        mlab.savefig("potential_xp_view.png")

        fig.scene.y_minus_view()
        mlab.savefig("potential_ym_view.png")

        fig.scene.y_plus_view()
        mlab.savefig("potential_yp_view.png")

        fig.scene.z_minus_view()
        mlab.savefig("potential_zm_view.png")

        fig.scene.z_plus_view()
        mlab.savefig("potential_zp_view.png")

    if interactive is True:
        # Enable interactive plot
        mlab.show()
    else:
        mlab.close(fig)




if __name__ == "__main__":
    iom = IOManager()

    # Read file with simulation data
    try:
        iom.open_file(filename=sys.argv[1])
    except IndexError:
        iom.open_file()

    parameters = iom.load_parameters()

    if parameters["dimension"] == 2:
        Potential = PotentialFactory().create_potential(parameters)
        Grid = TensorProductGrid(parameters)
        plot_potential(Grid, Potential, interactive=True)
    else:
        print("Not a potential in two space dimensions, silent return!")

    iom.finalize()
