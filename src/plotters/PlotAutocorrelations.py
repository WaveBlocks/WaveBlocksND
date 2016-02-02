"""The WaveBlocks Project

Plot the autocorrelations of the different wavepackets as well as the sum of all autocorrelations.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import max, real, imag, abs, add
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND.Plot import legend
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD
from functools import reduce


def read_all_datablocks(iom):
    """Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_autocorrelation(blockid=blockid):
            plot_autocorrelations(read_data(iom, blockid=blockid), blockid=blockid)
        else:
            print("Warning: Not plotting autocorrelations in block '%s'" % blockid)


def read_data(iom, blockid=0):
    """
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    if iom.has_parameters():
        parameters = iom.load_parameters()
        if "dt" in parameters:
            dt = parameters["dt"]
    else:
        dt = None

    timegrid = iom.load_autocorrelation_timegrid(blockid=blockid)

    autocorrelations = iom.load_autocorrelation(blockid=blockid, split=True)
    # Compute the sum over all levels
    autocorrelations.append(reduce(add, autocorrelations))

    return (timegrid, autocorrelations, dt)


def plot_autocorrelations(data, blockid=0, view=None):
    print("Plotting the autocorrelations of data block '%s'" % blockid)

    timegrid, autocorrelations, dt = data

    if dt is None:
        xlbl = r"Timesteps $n$"
        dt = 1.0
    else:
        xlbl = r"Time $t$"

    # View
    if view[0] is None:
        view[0] = dt * timegrid.min()
    if view[1] is None:
        view[1] = dt * timegrid.max()
    if view[2] is None:
        view[2] = 0.0
    if view[3] is None:
        view[3] = 1.1

    # Plot the autocorrelations
    fig = figure()
    ax = fig.gca()

    # Plot the autocorrelations of the individual wavepackets
    for i, datum in enumerate(autocorrelations[:-1]):
        ax.plot(timegrid*dt, abs(datum), label=r"$|\langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle|$" % (i,i))

    # Plot the sum of all autocorrelations
    ax.plot(timegrid*dt, abs(autocorrelations[-1]), color=(1,0,0), label=r"$\sum_i {|\langle \Phi_i(0) | \Phi_i(t) \rangle|}$")

    ax.grid(True)
    ax.set_xlim(view[:2])
    ax.set_ylim(view[2:])
    ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    ax.set_title(r"Autocorrelations of $\Psi$")
    legend(loc="upper right")
    ax.set_xlabel(xlbl)
    fig.savefig("autocorrelations_block"+str(blockid)+GD.output_format)
    close(fig)


    # Plot the autocorrelations
    N = len(autocorrelations) -1
    fig = figure()

    # Plot the autocorrelations of the individual wavepackets
    for i, datum in enumerate(autocorrelations[:-1]):
        ax = fig.add_subplot(N, 1, i+1)
        ax.plot(timegrid*dt, real(datum), label=r"$\Re \langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle$" % (i,i))
        ax.plot(timegrid*dt, imag(datum), label=r"$\Im \langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle$" % (i,i))
        ax.plot(timegrid*dt, abs(datum), label=r"$|\langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle|$" % (i,i))

        ax.grid(True)
        ax.set_xlim(view[:2])
        ax.set_ylim(view[2:])
        ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
        ax.set_xlabel(xlbl)
        legend(loc="upper right")

    ax.set_title(r"Autocorrelations of $\Psi$")
    fig.savefig("autocorrelations_per_component_block"+str(blockid)+GD.output_format)
    close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-b", "--blockid",
                        type = str,
                        help = "The data block to handle",
                        nargs = "*",
                        default = ["all"])

    parser.add_argument("-t", "--trange",
                        type = float,
                        help = "The plot range on the t-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-v", "--vrange",
                        type = float,
                        help = "The plot range on the y-axis",
                        nargs = 2,
                        default = [None, None])

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=args.datafile)

    # Which blocks to handle
    blockids = iom.get_block_ids()
    if not "all" in args.blockid:
        blockids = [ bid for bid in args.blockid if bid in blockids ]

    # The axes rectangle that is plotted
    view = args.trange + args.vrange

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting autocorrelations in data block '%s'" % blockid)

        if iom.has_autocorrelation(blockid=blockid):
            plot_autocorrelations(read_data(iom, blockid=blockid), blockid=blockid, view=view)
        else:
            print("Warning: Not plotting autocorrelations in block '%s'" % blockid)

    iom.finalize()
