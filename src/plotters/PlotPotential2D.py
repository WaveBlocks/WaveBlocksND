"""The WaveBlocks Project

Plot the eigenvalues (energy levels) of the potential.
This script is only for two-dimensional potentials.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import real

try:
    # New Mayavi version
    from mayavi import mlab
except ImportError:
    # Older Mayavi versions
    from enthought.mayavi import mlab

from WaveBlocksND import BlockFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager
from WaveBlocksND import ParameterLoader
from WaveBlocksND import GlobalDefaults as GLD


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
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-p", "--paramfile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = None)

    parser.add_argument("-b", "--blockid",
                        help = "The data block to handle",
                        nargs = "*",
                        default = [0])

    parser.add_argument("-x", "--xrange",
                        type = float,
                        help = "The plot range on the x-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-y", "--yrange",
                        type = float,
                        help = "The plot range on the y-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-z", "--zrange",
                        type = float,
                        help = "The plot range on the z-axis",
                        nargs = 2,
                        default = [0, 10])

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    try:
        iom.open_file(filename=args.datafile)
    except IndexError:
        iom.open_file()

    # Read file with parameter data for grid
    parameters = iom.load_parameters()

    if args.paramfile:
        PL = ParameterLoader()
        gridparams = PL.load_from_file(args.paramfile)
    else:
        gridparams = parameters

    # The axes rectangle that is plotted
    view = args.xrange + args.yrange + args.zrange
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
