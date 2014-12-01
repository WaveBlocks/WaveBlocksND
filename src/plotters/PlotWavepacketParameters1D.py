r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import real, imag, abs
from matplotlib.pyplot import figure, close

from WaveBlocksND import ComplexMath
from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD


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
    qhist, phist, Qhist, Phist, Shist = Pi

    # The Dimension D, we know that q has shape (#timesteps, D, 1)
    D = qhist.shape[1]

    if not D == 1:
        raise NotImplementedError("This script is for 1D wavepackets only")

    data = ( time,
             [qhist.reshape(-1)],
             [phist.reshape(-1)],
             [Qhist.reshape(-1)],
             [Phist.reshape(-1)],
             [Shist.reshape(-1)] )

    return data


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

    # The Dimension D, we know that q_0 has shape (#timesteps, D, 1)
    D = Pis[0][0].shape[1]
    if not D == 1:
        raise NotImplementedError("This script is for 1D wavepackets only")

    # List with Pi time evolutions
    Phist = []
    Qhist = []
    Shist = []
    phist = []
    qhist = []

    for q,p,Q,P,S in Pis:
        qhist.append(q.reshape(-1))
        phist.append(p.reshape(-1))
        Qhist.append(Q.reshape(-1))
        Phist.append(P.reshape(-1))
        Shist.append(S.reshape(-1))

    return (time, qhist, phist, Qhist, Phist, Shist)


def plot_parameters(data, index=0):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    print("Plotting the parameters of data block '%s'" % index)

    timegrid, qhist, phist, Qhist, Phist, Shist = data

    # Plot the time evolution of the parameters q, p, Q, P and S
    fig = figure(figsize=(12,12))

    ax = fig.add_subplot(4,2,1)
    for item in qhist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$q$")

    ax = fig.add_subplot(4,2,2)
    for item in phist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$p$")

    ax = fig.add_subplot(4,2,3)
    for item in Qhist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$\Re Q$")

    ax = fig.add_subplot(4,2,4)
    for item in Phist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$\Re P$")

    ax = fig.add_subplot(4,2,5)
    for item in Qhist:
        ax.plot(timegrid, imag(item))
    ax.grid(True)
    ax.set_title(r"$\Im Q$")

    ax = fig.add_subplot(4,2,6)
    for item in Phist:
        ax.plot(timegrid, imag(item))
    ax.grid(True)
    ax.set_title(r"$\Im P$")

    ax = fig.add_subplot(4,2,7)
    for item in Shist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig("wavepacket_parameters_block"+str(index)+GD.output_format)
    close(fig)


    # Plot the time evolution of the parameters P, Q, S, p and q
    # This time plot abs/angle instead of real/imag
    fig = figure(figsize=(12,12))

    ax = fig.add_subplot(4,2,1)
    for item in qhist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$q$")

    ax = fig.add_subplot(4,2,2)
    for item in phist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$p$")

    ax = fig.add_subplot(4,2,3)
    for item in Qhist:
        ax.plot(timegrid, abs(item))
    ax.grid(True)
    ax.set_title(r"$|Q|$")

    ax = fig.add_subplot(4,2,4)
    for item in Phist:
        ax.plot(timegrid, abs(item))
    ax.grid(True)
    ax.set_title(r"$|P|$")

    ax = fig.add_subplot(4,2,5)
    for item in Qhist:
        ax.plot(timegrid, ComplexMath.cont_angle(item))
    ax.grid(True)
    ax.set_title(r"$\arg Q$")

    ax = fig.add_subplot(4,2,6)
    for item in Phist:
        ax.plot(timegrid, ComplexMath.cont_angle(item))
    ax.grid(True)
    ax.set_title(r"$\arg P$")

    ax = fig.add_subplot(4,2,7)
    for item in Shist:
        ax.plot(timegrid, real(item))
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig("wavepacket_parameters_abs_ang_block"+str(index)+GD.output_format)
    close(fig)


    # Plot the complex trajectory of the parameters P
    fig = figure()
    ax = fig.gca()
    for item in Phist:
        ax.plot(real(item), imag(item), "-o")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_title(r"Trajectory of $P$")
    fig.savefig("wavepacket_parameters_trajectoryP_block"+str(index)+GD.output_format)
    close(fig)


    # Plot the complex trajectory of the parameters Q
    fig = figure()
    ax = fig.gca()
    for item in Qhist:
        ax.plot(real(item), imag(item), "-o")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_title(r"Trajectory of $Q$")
    fig.savefig("wavepacket_parameters_trajectoryQ_block"+str(index)+GD.output_format)
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
        print("Plotting wavepacket parameters in data block '%s'" % blockid)

        # NOTE: Add new algorithms here

        if iom.has_wavepacket(blockid=blockid):
            plot_parameters(read_data_homogeneous(iom, blockid=blockid), index=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_parameters(read_data_inhomogeneous(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting wavepacket parameters in block '%s'" % blockid)

    iom.finalize()
