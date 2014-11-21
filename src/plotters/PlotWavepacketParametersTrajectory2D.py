r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import squeeze
from matplotlib.pyplot import *

from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD


def read_all_datablocks(iom):
    r"""Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_wavepacket(blockid=blockid):
            plot_parameters(read_data_homogeneous(iom, blockid=blockid), index=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_parameters(read_data_inhomogeneous(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting wavepacket parameters in block '"+str(blockid)+"'!")


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if parameters.has_key("dt") else 1.0
    time = timegrid * dt

    Pi = iom.load_wavepacket_parameters(blockid=blockid)

    # The Dimension D, we know that q has shapw (#timesteps, D, 1)
    D = Pi[0].shape[1]
    if not D == 2:
        print("Warning: Trajectory plotting implemented only for 2D wavepackets")
        return

    Pi = map(squeeze, Pi)
    qhist, phist, Qhist, Phist, Shist = Pi
    return (time, [qhist], [phist], [Qhist], [Phist], [Shist])


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if parameters.has_key("dt") else 1.0
    time = timegrid * dt

    Pis = iom.load_inhomogwavepacket_parameters(blockid=blockid)

    # The Dimension D, we know that q_0 has shapw (#timesteps, D, 1)
    D = Pis[0][0].shape[1]
    if not D == 2:
        print("Warning: Trajectory plotting implemented only for 2D wavepackets")
        return

    # List with Pi time evolutions
    Phist = []
    Qhist = []
    Shist = []
    phist = []
    qhist = []

    for q,p,Q,P,S in Pis:
        qhist.append(squeeze(q))
        phist.append(squeeze(p))
        Qhist.append(squeeze(Q))
        Phist.append(squeeze(P))
        Shist.append(squeeze(S))

    return (time, qhist, phist, Qhist, Phist, Shist)


def plot_parameters(data, index=0):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    if data is None:
        return
    else:
        print("Plotting the parameters of data block '"+str(index)+"'")
        timegrid, qhist, phist, Qhist, Phist, Shist = data

    # Plot the 2D trajectory of the parameters q and p
    fig = figure()
    ax = fig.gca()
    for item in qhist:
        ax.plot(item[:,0], item[:,1], "-o", label=r"Trajectory of $q$")
    ax.set_xlabel(r"$q_x$")
    ax.set_ylabel(r"$q_y$")
    ax.grid(True)
    ax.set_title(r"Trajectory of $q$")
    fig.savefig("wavepacket_parameters_trajectoryq_block"+str(index)+GD.output_format)
    close(fig)

    fig = figure()
    ax = fig.gca()
    for item in phist:
        ax.plot(item[:,0], item[:,1], "-o", label=r"Trajectory of $p$")
    ax.set_xlabel(r"$p_x$")
    ax.set_ylabel(r"$p_y$")
    ax.grid(True)
    ax.set_title(r"Trajectory of $p$")
    fig.savefig("wavepacket_parameters_trajectoryp_block"+str(index)+GD.output_format)
    close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-b", "--blockid",
                        help = "The data block to handle",
                        nargs = "*",
                        default = [0])

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    try:
        iom.open_file(filename=args.datafile)
    except IndexError:
        iom.open_file()

    # Read the data and plot it, one plot for each data block.
    read_all_datablocks(iom)

    iom.finalize()
