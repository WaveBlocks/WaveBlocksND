#!/usr/bin/env python
r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real, imag, abs, where, nan, nanmin, nanmax, isnan, conjugate
from matplotlib.pyplot import figure, close

from WaveBlocksND import ComplexMath
from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Interface import GraphicsDefaults as GD


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if "dt" in parameters else 1.0
    # Filter
    time = timegrid * dt
    time = where(timegrid < 0, nan, time)

    Pi = iom.load_wavepacket_parameters(blockid=blockid)
    qhist, phist, Qhist, Phist, Shist = Pi

    # qhist, phist have shape (n, D, 1)
    # Qhist, Phist have shape (n, D, D)
    _, D, _ = Qhist.shape

    if not D == 1:
        raise NotImplementedError("This script is for 1D wavepackets only")

    return (time,
            [qhist.reshape(-1)],
            [phist.reshape(-1)],
            [Qhist.reshape(-1)],
            [Phist.reshape(-1)],
            [Shist.reshape(-1)])


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if "dt" in parameters else 1.0
    # Filter
    time = timegrid * dt
    time = where(timegrid < 0, nan, time)

    Pis = iom.load_inhomogwavepacket_parameters(blockid=blockid)

    # qhist, phist have shape (n, D, 1)
    # Qhist, Phist have shape (n, D, D)
    _, D, _ = Pis[0][0].shape

    if not D == 1:
        raise NotImplementedError("This script is for 1D wavepackets only")

    # List with Pi time evolutions
    Phist = []
    Qhist = []
    Shist = []
    phist = []
    qhist = []

    for q, p, Q, P, S in Pis:
        qhist.append(q.reshape(-1))
        phist.append(p.reshape(-1))
        Qhist.append(Q.reshape(-1))
        Phist.append(P.reshape(-1))
        Shist.append(S.reshape(-1))

    return (time, qhist, phist, Qhist, Phist, Shist)


def plot_parameters(data, index=0, view=[None, None], path='.'):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    print("Plotting the parameters of data block '{}'".format(index))

    time, qhist, phist, Qhist, Phist, Shist = data

    # View
    if view[0] is None:
        view[0] = nanmin(time)
    if view[1] is None:
        view[1] = nanmax(time)

    # Plot the time evolution of the parameters q, p, Q, P and S
    fig = figure(figsize=(12, 12))

    ax = fig.add_subplot(4, 2, 1)
    for item in qhist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$q$")

    ax = fig.add_subplot(4, 2, 2)
    for item in phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$p$")

    ax = fig.add_subplot(4, 2, 3)
    for item in Qhist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Re Q$")

    ax = fig.add_subplot(4, 2, 4)
    for item in Phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Re P$")

    ax = fig.add_subplot(4, 2, 5)
    for item in Qhist:
        ax.plot(time, imag(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Im Q$")

    ax = fig.add_subplot(4, 2, 6)
    for item in Phist:
        ax.plot(time, imag(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Im P$")

    ax = fig.add_subplot(4, 2, 7)
    for item in Shist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig(os.path.join(path, "wavepacket_parameters_block".format(index, GD.output_format)))
    close(fig)


    # Plot the time evolution of the parameters P, Q, S, p and q
    # This time plot abs/angle instead of real/imag
    fig = figure(figsize=(12, 12))

    ax = fig.add_subplot(4, 2, 1)
    for item in qhist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$q$")

    ax = fig.add_subplot(4, 2, 2)
    for item in phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$p$")

    ax = fig.add_subplot(4, 2, 3)
    for item in Qhist:
        ax.plot(time, abs(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$|Q|$")

    ax = fig.add_subplot(4, 2, 4)
    for item in Phist:
        ax.plot(time, abs(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$|P|$")

    ax = fig.add_subplot(4, 2, 5)
    for item in Qhist:
        ax.plot(time, ComplexMath.cont_angle(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\arg Q$")

    ax = fig.add_subplot(4, 2, 6)
    for item in Phist:
        ax.plot(time, ComplexMath.cont_angle(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\arg P$")

    ax = fig.add_subplot(4, 2, 7)
    for item in Shist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig(os.path.join(path, "wavepacket_parameters_abs_ang_block{}{}".format(index, GD.output_format)))
    close(fig)


    vals = ~isnan(time)

    # Plot the complex trajectory of the parameters P
    fig = figure()
    ax = fig.gca()
    for item in Phist:
        ax.plot(real(item[vals]), imag(item[vals]), "-o")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Re P(t)$")
    ax.set_ylabel(r"$\Im P(t)$")
    ax.set_title(r"Trajectory of $P$")
    fig.savefig(os.path.join(path, "wavepacket_parameters_trajectoryP_block{}{}".format(index, GD.output_format)))
    close(fig)


    # Plot the complex trajectory of the parameters Q
    fig = figure()
    ax = fig.gca()
    for item in Qhist:
        ax.plot(real(item[vals]), imag(item[vals]), "-o")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Re Q(t)$")
    ax.set_ylabel(r"$\Im Q(t)$")
    ax.set_title(r"Trajectory of $Q$")
    fig.savefig(os.path.join(path, "wavepacket_parameters_trajectoryQ_block{}{}".format(index, GD.output_format)))
    close(fig)


    fig = figure()
    ax = fig.gca()
    for P, Q in zip(Phist, Qhist):
        ax.plot(time, abs(P.T * Q - Q.T * P), label=r'$P^{T} Q - Q^{T} P = 0$')
        ax.plot(time, abs(conjugate(P.T) * Q - conjugate(Q.T) * P + 2.0j), label=r'$P^{H} Q - Q^{H} P + 2\imath I = 0$')
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_xlabel(r"$t$")
    ax.legend(loc='upper left')
    fig.suptitle("Wavepacket PQ relation check")
    fig.savefig(os.path.join(path, "wavepacket_PQrelation_block{}{}".format(index, GD.output_format)))
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

    parser.add_argument("-r", "--resultspath",
                        type = str,
                        help = "Path where to put the results.",
                        nargs = "?",
                        default = '.')

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

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting wavepacket parameters in data block '{}'".format(blockid))

        # NOTE: Add new algorithms here

        if iom.has_wavepacket(blockid=blockid):
            plot_parameters(read_data_homogeneous(iom, blockid=blockid), index=blockid, path=resultspath)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_parameters(read_data_inhomogeneous(iom, blockid=blockid), index=blockid, path=resultspath)
        else:
            print("Warning: Not plotting wavepacket parameters in block '{}'".format(blockid))

    iom.finalize()
