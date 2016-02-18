#!/usr/bin/env python
"""The WaveBlocks Project

Plot the eigenvalues (energy levels) of the potential.
This script is only for one-dimensional potentials.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real
from matplotlib.pyplot import figure, close

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager
from WaveBlocksND.Plot import legend
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Interface import GraphicsDefaults as GD


def plot_potential(grid, potential, view=None, size=(12, 9), path='.'):
    # The Grid
    u = grid.get_nodes(split=True)
    u = real(u[0])

    # Create potential and evaluate eigenvalues
    potew = potential.evaluate_eigenvalues_at(grid)
    potew = [real(level).reshape(-1) for level in potew]

    # View
    if view[0] is None:
        view[0] = u.min()
    if view[1] is None:
        view[1] = u.max()

    # Plot the energy surfaces of the potential
    fig = figure(figsize=size)
    ax = fig.gca()

    for index, ew in enumerate(potew):
        ax.plot(u, ew, label=r"$\lambda_{%d}$" % index)

    ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.grid(True)
    ax.set_xlim(view[:2])
    ax.set_ylim(view[2:])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\lambda_i\left(x\right)$")
    legend(loc="outer right")
    ax.set_title(r"The eigenvalues $\lambda_i$ of the potential $V\left(x\right)$")
    fig.savefig(os.path.join(path, "potential"+GD.output_format))
    close(fig)




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
                        default = [-2, 5])

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
    view = args.xrange + args.yrange

    if parameters["dimension"] == 1:
        Potential = BlockFactory().create_potential(parameters)
        Grid = TensorProductGrid(gridparams["limits"], gridparams["number_nodes"])
        plot_potential(Grid, Potential, view=view, path=resultspath)
    else:
        print("Not a potential in one space dimension, silent return!")

    iom.finalize()
