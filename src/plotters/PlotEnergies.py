"""The WaveBlocks Project

Plot the energies of the different wavepackets as well as the sum of all energies.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import abs, add
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND.Plot import legend
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD


def read_all_datablocks(iom):
    """Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """


def read_data(iom, blockid=0):
    """
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    timegridk, timegridp = iom.load_energy_timegrid(blockid=blockid)
    have_dt = iom.has_parameters()
    if have_dt:
        parameters = iom.load_parameters()
        timek = timegridk * parameters["dt"]
        timep = timegridp * parameters["dt"]
    else:
        timek = timegridk
        timep = timegridp

    ekin, epot = iom.load_energy(blockid=blockid, split=True)

    # Compute the sum of all energies
    ekin.append(reduce(add, ekin))
    epot.append(reduce(add, epot))

    return (timek, timep, ekin, epot, have_dt)


def plot_energies(data, index=0):
    print("Plotting the energies of data block '%s'" % index)

    timegridk, timegridp, ekin, epot, have_dt = data

    xlbl = r"Time $t$" if have_dt else r"Timesteps $n$"

    # Plot the energies
    fig = figure()
    ax = fig.gca()

    # Plot the kinetic energy of the individual wave packets
    for i, kin in enumerate(ekin[:-1]):
        ax.plot(timegridk, kin, label=r"$E^{kin}_{%d}$" % i)

    # Plot the potential energy of the individual wave packets
    for i, pot in enumerate(epot[:-1]):
        ax.plot(timegridp, pot, label=r"$E^{pot}_{%d}$" % i)

    # Plot the sum of kinetic and potential energy for all wave packets
    for i, (kin, pot) in enumerate(zip(ekin, epot)[:-1]):
        ax.plot(timegridk, kin + pot, label=r"$E^{kin}_{%d}+E^{pot}_{%d}$" % (i,i))

    # Plot sum of kinetic and sum of potential energy
    ax.plot(timegridk, ekin[-1], label=r"$\sum_i E^{kin}_i$")
    ax.plot(timegridp, epot[-1], label=r"$\sum_i E^{pot}_i$")

    # Plot the overall energy of all wave packets
    ax.plot(timegridk, ekin[-1] + epot[-1], label=r"$\sum_i E^{kin}_i + \sum_i E^{pot}_i$")

    ax.set_xlim(min(min(timegridk),min(timegridp)), max(max(timegridk),max(timegridp)))
    ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    ax.grid(True)
    ax.set_xlabel(xlbl)
    legend(loc="outer right")
    ax.set_title(r"Energies of the wavepacket $\Psi$")

    fig.savefig("energies_block"+str(index)+GD.output_format)
    close(fig)


    # Plot the energy drift
    e_orig = (ekin[-1]+epot[-1])[0]
    data = abs(e_orig-(ekin[-1]+epot[-1]))

    fig = figure()
    ax = fig.gca()

    ax.plot(timegridk, data, label=r"$|E_O^0 - \left( E_k^0 + E_p^0 \right) |$")

    ax.set_xlim(min(timegridk), max(timegridk))
    ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
    ax.grid(True)
    ax.set_xlabel(xlbl)
    ax.set_ylabel(r"$|E_O^0 - \left( E_k^0 + E_p^0 \right) |$")
    ax.set_title(r"Energy drift of the wavepacket $\Psi$")

    fig.savefig("energy_drift_block"+str(index)+"_lin"+GD.output_format)
    close(fig)


    fig = figure()
    ax = fig.gca()

    ax.semilogy(timegridk, data, label=r"$|E_O^0 - \left( E_k^0 + E_p^0 \right) |$")

    ax.set_xlim(min(timegridk), max(timegridk))
    ax.grid(True)
    ax.set_xlabel(xlbl)
    ax.set_ylabel(r"$|E_O^0 - \left( E_k^0 + E_p^0 \right) |$")
    ax.set_title(r"Energy drift of the wavepacket $\Psi$")

    fig.savefig("energy_drift_block"+str(index)+"_log"+GD.output_format)
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
    for blockid in iom.get_block_ids():
        print("Plotting energies in data block '%s'" % blockid)

        if iom.has_energy(blockid=blockid):
            plot_energies(read_data(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting energies in block '%s'" % blockid)

    iom.finalize()
