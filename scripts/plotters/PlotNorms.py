#!/usr/bin/env python
"""The WaveBlocks Project

Plot the norms of the different wavepackets as well as the sum of all norms.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from functools import reduce
from numpy import add, square, sqrt, max, where, nan, nanmin, nanmax
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND.Plot import legend
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Interface import GraphicsDefaults as GD


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

    timegrid = iom.load_norm_timegrid(blockid=blockid)

    norms = iom.load_norm(blockid=blockid, split=True)
    # Compute the sum over all levels
    norms.append(sqrt(reduce(add, list(map(square, norms)))))

    return (timegrid, norms, dt)


def plot_norms(data, blockid=0, view=None, path='.'):
    print("Plotting the norms of data block '{}'".format(blockid))

    timegrid, norms, dt = data

    if dt is None:
        xlbl = r"Timesteps $n$"
        dt = 1.0
    else:
        xlbl = r"Time $t$"

    # Filter
    time = timegrid * dt
    time = where(timegrid < 0, nan, time)

    # View
    if view[0] is None:
        view[0] = nanmin(time)
    if view[1] is None:
        view[1] = nanmax(time)
    if view[2] is None:
        view[2] = 0.0
    if view[3] is None:
        view[3] = 1.1 * max(norms[-1])

    # Plot the norms
    fig = figure()
    ax = fig.gca()

    # Plot the norms of the individual wavepackets
    for i, norm in enumerate(norms[:-1]):
        label_i = r"$\| \Phi_{%d} \|$" % i
        ax.plot(time, norm, label=label_i)

    # Plot the sum of all norms
    ax.plot(time, norms[-1], color=(1, 0, 0), label=r"${\sqrt{\sum_i {\| \Phi_i \|^2}}}$")

    ax.grid(True)
    ax.set_xlim(view[:2])
    ax.set_ylim(view[2], view[3])
    ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.set_title(r"Norms of $\Psi$")
    legend(loc="outer right")
    ax.set_xlabel(xlbl)
    fig.savefig(os.path.join(path, "norms_block" + str(blockid) + GD.output_format))
    close(fig)


    fig = figure()
    ax = fig.gca()

    # Plot the squared norms of the individual wavepackets
    for i, norm in enumerate(norms[:-1]):
        label_i = r"$\| \Phi_{%d} \|^2$" % i
        ax.plot(time, norm**2, label=label_i)

    # Plot the squared sum of all norms
    ax.plot(time, norms[-1]**2, color=(1, 0, 0), label=r"${\sum_i {\| \Phi_i \|^2}}$")

    ax.grid(True)
    ax.set_xlim(view[:2])
    ax.set_ylim(view[2], view[3]**2)
    ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.set_title(r"Squared norms of $\Psi$")
    legend(loc="outer right")
    ax.set_xlabel(xlbl)
    fig.savefig(os.path.join(path, "norms_sqr_block" + str(blockid) + GD.output_format))
    close(fig)


    # Plot the difference from the theoretical norm
    fig = figure()
    ax = fig.gca()

    ax.plot(time, abs(norms[-1][0] - norms[-1]), label=r"$\|\Psi\|_0 - \|\Psi\|_t$")

    ax.grid(True)
    ax.set_xlim(view[:2])
    ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.set_title(r"Drift of $\| \Psi \|$")
    legend(loc="outer right")
    ax.set_xlabel(xlbl)
    ax.set_ylabel(r"$\|\Psi\|_0 - \|\Psi\|_t$")
    fig.savefig(os.path.join(path, "norms_drift_block" + str(blockid) + GD.output_format))
    close(fig)


    fig = figure()
    ax = fig.gca()

    ax.semilogy(time, abs(norms[-1][0] - norms[-1]), label=r"$\|\Psi\|_0 - \|\Psi\|_t$")

    ax.set_xlim(view[:2])
    ax.grid(True)
    ax.set_title(r"Drift of $\| \Psi \|$")
    legend(loc="outer right")
    ax.set_xlabel(xlbl)
    ax.set_ylabel(r"$\|\Psi\|_0 - \|\Psi\|_t$")
    fig.savefig(os.path.join(path, "norms_drift_block" + str(blockid) + "_log" + GD.output_format))
    close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type=str,
                        help="The simulation data file",
                        nargs="?",
                        default=GLD.file_resultdatafile)

    parser.add_argument("-b", "--blockid",
                        type=str,
                        help="The data block to handle",
                        nargs="*",
                        default=["all"])

    parser.add_argument("-r", "--resultspath",
                        type=str,
                        help="Path where to put the results.",
                        nargs="?",
                        default='.')

    parser.add_argument("-t", "--trange",
                        type=float,
                        help="The plot range on the t-axis",
                        nargs=2,
                        default=[None, None])

    parser.add_argument("-v", "--vrange",
                        type=float,
                        help="The plot range on the y-axis",
                        nargs=2,
                        default=[None, None])

    args = parser.parse_args()


    # File with the simulation data
    resultspath = os.path.abspath(args.resultspath)

    if not os.path.exists(resultspath):
        raise IOError("The results path does not exist: {}".format(args.resultspath))

    datafile = os.path.abspath(os.path.join(args.resultspath, args.datafile))

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=datafile)

    # Which blocks to handle
    blockids = iom.get_block_ids()
    if "all" not in args.blockid:
        blockids = [bid for bid in args.blockid if bid in blockids]

    # The axes rectangle that is plotted
    view = args.trange + args.vrange

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting norms in data block '{}'".format(blockid))

        if iom.has_norm(blockid=blockid):
            plot_norms(read_data(iom, blockid=blockid), blockid=blockid, view=view, path=resultspath)
        else:
            print("Warning: Not plotting norms in block '{}'".format(blockid))

    iom.finalize()
