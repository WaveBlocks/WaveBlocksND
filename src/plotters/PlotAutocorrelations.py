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


def read_all_datablocks(iom):
    """Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_autocorrelation(blockid=blockid):
            plot_autocorrelations(read_data(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting autocorrelations in block '%s'" % blockid)


def read_data(iom, blockid=0):
    """
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    have_dt = iom.has_parameters()
    if have_dt:
        parameters = iom.load_parameters()

    dt = parameters["dt"] if parameters.has_key("dt") else 1.0
    timegrid = iom.load_autocorrelation_timegrid(blockid=blockid)
    time = timegrid * dt

    autocorrelations = iom.load_autocorrelation(blockid=blockid, split=True)
    autocorrelations.append(reduce(add, autocorrelations))

    return (time, autocorrelations, have_dt)


def plot_autocorrelations(data, index=0):
    print("Plotting the autocorrelations of data block '%s'" % index)

    timegrid, autocorrelations, have_dt = data

    xlbl = r"Time $t$" if have_dt else r"Timesteps $n$"

    # Plot the autocorrelations
    fig = figure()
    ax = fig.gca()

    # Plot the autocorrelations of the individual wavepackets
    for i, datum in enumerate(autocorrelations[:-1]):
        ax.plot(timegrid, abs(datum), label=r"$|\langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle|$" % (i,i))

    # Plot the sum of all autocorrelations
    ax.plot(timegrid, abs(autocorrelations[-1]), color=(1,0,0), label=r"$\sum_i {|\langle \Phi_i(0) | \Phi_i(t) \rangle|}$")

    ax.set_xlim(min(timegrid), max(timegrid))
    ax.grid(True)
    ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    ax.set_title(r"Autocorrelations of $\Psi$")
    legend(loc="upper right")
    ax.set_xlabel(xlbl)
    ax.set_ylim(0, 1.1)
    fig.savefig("autocorrelations_block"+str(index)+GD.output_format)
    close(fig)


    # Plot the autocorrelations
    N = len(autocorrelations) -1
    fig = figure()

    # Plot the autocorrelations of the individual wavepackets
    for i, datum in enumerate(autocorrelations[:-1]):
        ax = fig.add_subplot(N, 1, i+1)
        ax.plot(timegrid, real(datum), label=r"$\Re \langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle$" % (i,i))
        ax.plot(timegrid, imag(datum), label=r"$\Im \langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle$" % (i,i))
        ax.plot(timegrid, abs(datum), label=r"$|\langle \Phi_{%d}(0) | \Phi_{%d}(t) \rangle|$" % (i,i))

        ax.set_xlim(min(timegrid), max(timegrid))
        ax.grid(True)
        ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
        ax.set_xlabel(xlbl)

        legend(loc="upper right")
    ax.set_ylim(0, 1.1)
    ax.set_title(r"Autocorrelations of $\Psi$")
    fig.savefig("autocorrelations_per_component_block"+str(index)+GD.output_format)
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

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=args.datafile)

    # Which blocks to handle
    if "all" in args.blockid:
        blockids = iom.get_block_ids()
    else:
        blockids = args.blockid

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting autocorrelations in data block '%s'" % blockid)

        if iom.has_autocorrelation(blockid=blockid):
            plot_autocorrelations(read_data(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting autocorrelations in block '%s'" % blockid)

    iom.finalize()
