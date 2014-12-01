r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import array, linspace, dot, sin, cos, pi, real
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    Pi = iom.load_wavepacket_parameters(blockid=blockid)
    qhist, phist, Qhist, Phist, Shist = Pi

    # The Dimension D, we know that q has shape (#timesteps, D, 1)
    D = qhist.shape[1]
    if not D == 2:
        raise NotImplementedError("Trajectory plotting implemented only for 2D wavepackets")

    Pi = [ [real(qhist.reshape((-1,D)))],
           [real(phist.reshape((-1,D)))],
           [Qhist.reshape((-1,D,D))],
           [Phist.reshape((-1,D,D))],
           [Shist.reshape((-1,1))] ]

    return Pi


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    Pis = iom.load_inhomogwavepacket_parameters(blockid=blockid)

    # The Dimension D, we know that q_0 has shape (#timesteps, D, 1)
    D = Pis[0][0].shape[1]
    if not D == 2:
        raise NotImplementedError("Trajectory plotting implemented only for 2D wavepackets")

    # List with Pi time evolutions
    Phist = []
    Qhist = []
    Shist = []
    phist = []
    qhist = []

    for q,p,Q,P,S in Pis:
        qhist.append(real(q.reshape((-1,D,))))
        phist.append(real(p.reshape((-1,D,))))
        Qhist.append(Q.reshape((-1,D,D)))
        Phist.append(P.reshape((-1,D,D)))
        Shist.append(S.reshape((-1,1,)))

    return (qhist, phist, Qhist, Phist, Shist)


def plot_parameters(data, index=0):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    print("Plotting the parameters of data block '%s'" % index)

    qhist, phist, Qhist, Phist, Shist = data

    r = 1.0
    theta = linspace(0, 2*pi, 64)
    circle = array([r*cos(theta), r*sin(theta)])

    # Plot the 2D trajectory of the parameters q and p
    fig = figure()
    ax = fig.gca()
    for qitem, Qitem in zip(qhist, Qhist):
        for q, Q in zip(qitem, Qitem):
            C = real(q.reshape((-1,1)) + dot(Q, circle))
            ax.plot(q[0], q[1], "bo")
            ax.plot(C[0,:], C[1,:], "b-", alpha=0.5)

    ax.set_xlabel(r"$q_x$")
    ax.set_ylabel(r"$q_y$")
    ax.grid(True)
    ax.set_title(r"Schematic evolution of $q$ and $Q$")
    fig.savefig("wavepacket_parameters_schema_pos_block"+str(index)+GD.output_format)
    close(fig)


    fig = figure()
    ax = fig.gca()
    for pitem, Pitem in zip(phist, Phist):
        for p, P in zip(pitem, Pitem):
            C = real(p.reshape((-1,1)) + dot(P, circle))
            ax.plot(p[0], p[1], "bo")
            ax.plot(C[0,:], C[1,:], "b-", alpha=0.5)
    ax.set_xlabel(r"$p_x$")
    ax.set_ylabel(r"$p_y$")
    ax.grid(True)
    ax.set_title(r"Schematic evolution of $p$ and $P$")
    fig.savefig("wavepacket_parameters_schema_mom_lock"+str(index)+GD.output_format)
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
