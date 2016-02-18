#!/usr/bin/env python
"""The WaveBlocks Project

Plot the eigenvalues (energy levels) of the potential.
This script is only for two-dimensional potentials.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real
from mayavi import mlab

from WaveBlocksND import BlockFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager
from WaveBlocksND import ParameterLoader
from WaveBlocksND import GlobalDefaults as GLD


def plot_potential(grid, potential, sparsify=1, along_axes=False, view=None, interactive=False, path='.'):
    """Plot the potential

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # The Grid
    u, v = grid.get_nodes(split=True, flat=False)
    u = real(u[::sparsify, ::sparsify])
    v = real(v[::sparsify, ::sparsify])

    # Create potential and evaluate eigenvalues
    potew = potential.evaluate_eigenvalues_at(grid)
    potew = [real(level).reshape(grid.get_number_nodes(overall=False))[::sparsify, ::sparsify] for level in potew]

    # Plot
    if not interactive:
        mlab.options.offscreen = True

    fig = mlab.figure(size=(800, 700))

    for level in potew:
        # The energy surfaces of the potential
        src = mlab.pipeline.grid_source(u, v, level)

        # Clip to given view
        if view is not None:
            geometry_filter = mlab.pipeline.user_defined(src, filter='GeometryFilter')
            geometry_filter.filter.extent_clipping = True
            geometry_filter.filter.extent = view
            src = mlab.pipeline.user_defined(geometry_filter, filter='CleanPolyData')

        # Plot the surface
        normals = mlab.pipeline.poly_data_normals(src)
        mlab.pipeline.surface(normals)

    mlab.axes()

    fig.scene.parallel_projection = True
    fig.scene.isometric_view()
    # fig.scene.show_axes = True

    mlab.draw()
    if interactive:
        mlab.show()
    else:
        mlab.savefig(os.path.join(path, "potential_3D_view.png"))
        mlab.close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-p", "--parametersfile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = None)

    # parser.add_argument("-b", "--blockid",
    #                     type = str,
    #                     help = "The data block to handle",
    #                     nargs = "*",
    #                     default = ["all"])

    parser.add_argument("-r", "--resultspath",
                        type = str,
                        help = "Path where to put the results.",
                        nargs = "?",
                        default = '.')

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

    parser.add_argument("-s", "--sparsify",
                        type = int,
                        help = "Plot only every s-th point",
                        default = 10)

    parser.add_argument("-i", "--interactive",
                        action = "store_true",
                        help = "Hold the plot open for interactive manipulation")

    args = parser.parse_args()


    # File with the simulation data
    resultspath = os.path.abspath(args.resultspath)

    if not os.path.exists(resultspath):
        raise IOError("The results path does not exist: {}".format(args.resultspath))

    datafile = os.path.abspath(os.path.join(args.resultspath, args.datafile))
    parametersfile = os.path.abspath(os.path.join(args.resultspath, args.parametersfile))

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=datafile)

    # Read file with parameter data for grid
    parameters = iom.load_parameters()

    if args.parametersfile:
        PL = ParameterLoader()
        gridparams = PL.load_from_file(parametersfile)
    else:
        gridparams = parameters

    # The axes rectangle that is plotted
    view = args.xrange + args.yrange + args.zrange

    # Plot
    if parameters["dimension"] == 2:
        Potential = BlockFactory().create_potential(parameters)
        Grid = TensorProductGrid(gridparams["limits"], gridparams["number_nodes"])
        plot_potential(Grid, Potential,
                       sparsify=args.sparsify,
                       view=view,
                       interactive=args.interactive,
                       path=resultspath)
    else:
        print("Not a potential in two space dimensions, silent return!")

    iom.finalize()
