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

from WaveBlocksND import BlockFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager


def plot_potential(grid, potential, along_axes=False, interactive=False, view=None, size=(800,700)):
    # The Grid
    u, v = grid.get_nodes(split=True, flat=False)
    u = real(u)
    v = real(v)

    # Create potential and evaluate eigenvalues
    potew = potential.evaluate_eigenvalues_at(grid)
    potew = [ level.reshape(grid.get_number_nodes(overall=False)) for level in potew ]

    # Plot the energy surfaces of the potential
    fig = mlab.figure(size=size)

    for level in potew:
        if view is not None:
            mlab.surf(u, v, real(level), extent=view)
        else:
            mlab.surf(u, v, real(level))

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

    # Read file with parameter data for grid
    try:
        gridparams = PL.load_from_file(sys.argv[2])
    except IndexError:
        gridparams = parameters

    # Manually adjust the plotting region
    # view = [xmin, xmax, ymin, ymax, zmin, zmax]
    view = [None, None, None, None, 0, 10]

    if view[0] is None:
        view[0] = gridparams["limits"][0][0]
    if view[1] is None:
        view[1] = gridparams["limits"][0][1]
    if view[2] is None:
        view[2] = gridparams["limits"][1][0]
    if view[3] is None:
        view[3] = gridparams["limits"][1][1]

    print("Plotting in region: "+str(view))

    if parameters["dimension"] == 2:
        Potential = BlockFactory().create_potential(parameters)
        Grid = TensorProductGrid(gridparams["limits"], gridparams["number_nodes"])
        plot_potential(Grid, Potential, interactive=True, view=view)
    else:
        print("Not a potential in two space dimensions, silent return!")

    iom.finalize()
