#!/usr/bin/env python
r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real, imag, abs, array, where, nan, isnan, nanmin, nanmax
from numpy.linalg import norm, det
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

    qhist = array([norm(qhist[i, :]) for i in range(qhist.shape[0])]).reshape(-1)
    phist = array([norm(phist[i, :]) for i in range(phist.shape[0])]).reshape(-1)
    Qhist = array([det(Qhist[i, :, :]) for i in range(Qhist.shape[0])]).reshape(-1)
    Phist = array([det(Phist[i, :, :]) for i in range(Phist.shape[0])]).reshape(-1)
    Shist = Shist.reshape(-1)

    return (time, [qhist], [phist], [Qhist], [Phist], [Shist])


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

    # List with Pi time evolutions
    Phist = []
    Qhist = []
    Shist = []
    phist = []
    qhist = []

    for q, p, Q, P, S in Pis:
        qhist.append(array([norm(q[i, :]) for i in range(q.shape[0])]).reshape(-1))
        phist.append(array([norm(p[i, :]) for i in range(p.shape[0])]).reshape(-1))
        Qhist.append(array([det(Q[i, :, :]) for i in range(Q.shape[0])]).reshape(-1))
        Phist.append(array([det(P[i, :, :]) for i in range(P.shape[0])]).reshape(-1))
        Shist.append(S.reshape(-1))

    return (time, qhist, phist, Qhist, Phist, Shist)


def plot_parameters(data, index=0, view=[None, None], path='.'):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    print("Plotting the parameters of data block '%s'" % index)

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
    ax.set_title(r"$\|q\|_2$")

    ax = fig.add_subplot(4, 2, 2)
    for item in phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\|p\|_2$")

    ax = fig.add_subplot(4, 2, 3)
    for item in Qhist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Re \det Q$")

    ax = fig.add_subplot(4, 2, 4)
    for item in Phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Re \det P$")

    ax = fig.add_subplot(4, 2, 5)
    for item in Qhist:
        ax.plot(time, imag(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Im \det Q$")

    ax = fig.add_subplot(4, 2, 6)
    for item in Phist:
        ax.plot(time, imag(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\Im \det P$")

    ax = fig.add_subplot(4, 2, 7)
    for item in Shist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig(os.path.join(path, "wavepacket_parameters_block"+str(index)+GD.output_format))
    close(fig)


    # Plot the time evolution of the parameters P, Q, S, p and q
    # This time plot abs/angle instead of real/imag
    fig = figure(figsize=(12, 12))

    ax = fig.add_subplot(4, 2, 1)
    for item in qhist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\|q\|_2$")

    ax = fig.add_subplot(4, 2, 2)
    for item in phist:
        ax.plot(time, real(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\|p\|_2$")

    ax = fig.add_subplot(4, 2, 3)
    for item in Qhist:
        ax.plot(time, abs(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$|\det Q|$")

    ax = fig.add_subplot(4, 2, 4)
    for item in Phist:
        ax.plot(time, abs(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$|\det P|$")

    ax = fig.add_subplot(4, 2, 5)
    for item in Qhist:
        ax.plot(time, ComplexMath.cont_angle(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\arg \det Q$")

    ax = fig.add_subplot(4, 2, 6)
    for item in Phist:
        ax.plot(time, ComplexMath.cont_angle(item))
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$\arg \det P$")

    ax = fig.add_subplot(4, 2, 7)
    for item in Shist:
        ax.plot(time, real(item), label=r"$S$")
    ax.set_xlim(view[0], view[1])
    ax.grid(True)
    ax.set_title(r"$S$")

    fig.suptitle("Wavepacket parameters")
    fig.savefig(os.path.join(path, "wavepacket_parameters_abs_ang_block"+str(index)+GD.output_format))
    close(fig)


    vals = ~isnan(time)

    # Plot the complex trajectory of the parameters P
    fig = figure()
    ax = fig.gca()
    for item in Phist:
        ax.plot(real(item[vals]), imag(item[vals]), "-o")
    ax.set_xlabel(r"$\Re \det P(t)$")
    ax.set_ylabel(r"$\Im \det P(t)$")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_title(r"Trajectory of $\det P$")
    fig.savefig(os.path.join(path, "wavepacket_parameters_trajectoryP_block"+str(index)+GD.output_format))
    close(fig)


    # Plot the complex trajectory of the parameters Q
    fig = figure()
    ax = fig.gca()
    for item in Qhist:
        ax.plot(real(item[vals]), imag(item[vals]), "-o")
    ax.set_xlabel(r"$\Re \det Q(t)$")
    ax.set_ylabel(r"$\Im \det Q(t)$")
    ax.grid(True)
    ax.set_aspect("equal")
    ax.set_title(r"Trajectory of $\det Q$")
    fig.savefig(os.path.join(path, "wavepacket_parameters_trajectoryQ_block"+str(index)+GD.output_format))
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
        raise IOError("The results path does not exist: " + args.resultspath)

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
        print("Plotting wavepacket parameters in data block '%s'" % blockid)

        # NOTE: Add new algorithms here

        if iom.has_wavepacket(blockid=blockid):
            plot_parameters(read_data_homogeneous(iom, blockid=blockid), index=blockid, path=resultspath)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_parameters(read_data_inhomogeneous(iom, blockid=blockid), index=blockid, path=resultspath)
        else:
            print("Warning: Not plotting wavepacket parameters in block '%s'" % blockid)

    iom.finalize()
