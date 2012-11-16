"""The WaveBlocks Project

Plot the eigenvalues (energy levels) of the potential.
This script is only for one-dimensional potentials.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from numpy import real
from matplotlib.pyplot import *

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import TensorProductGrid
from WaveBlocksND import IOManager

import GraphicsDefaults as GD


def plot_potential(grid, potential, view=None, size=(12,9), fill=False):
    # The Grid
    u = grid.get_nodes(split=True, flat=False)
    u = real(u[0])

    # Create potential and evaluate eigenvalues
    potew = potential.evaluate_eigenvalues_at(grid)
    potew = [ level.reshape(grid.get_number_nodes(overall=False)) for level in potew ]

    # Plot the energy surfaces of the potential
    fig = figure(figsize=size)
    ax = fig.gca()

    for index, ew in enumerate(potew):
        if fill:
            ax.fill(u, ew, facecolor="blue", alpha=0.25)
        ax.plot(u, ew, label=r"$\lambda_"+str(index)+r"$")

    ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    ax.grid(True)
    # Set the aspect window
    if view is not None:
        ax.set_xlim(view[:2])
        ax.set_ylim(view[2:])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\lambda_i\left(x\right)$")
    legend(loc="outer right")
    ax.set_title(r"The eigenvalues $\lambda_i$ of the potential $V\left(x\right)$")
    fig.savefig("potential"+GD.output_format)
    close(fig)




if __name__ == "__main__":
    iom = IOManager()
    PL = ParameterLoader()

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
    # view = [xmin, xmax, ymin, ymax]
    view = [None, None, -2, 5]

    if view[0] is None:
        view[0] = gridparams["limits"][0][0]
    if view[1] is None:
        view[1] = gridparams["limits"][0][1]

    print("Plotting in region: "+str(view))

    if parameters["dimension"] == 1:
        Potential = BlockFactory().create_potential(parameters)
        Grid = TensorProductGrid(gridparams["limits"], gridparams["number_nodes"])
        plot_potential(Grid, Potential, view)
    else:
        print("Not a potential in one space dimension, silent return!")

    iom.finalize()
