#!/usr/bin/env python
"""The WaveBlocks Project

Plot the coefficients $c_i$ of each component of a homogeneous or
inhomogeneous Hagedorn wavepacket for all timesteps during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import abs, angle, array
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory
from WaveBlocksND.Plot import stemcf
from WaveBlocksND import GlobalDefaults as GLD


def read_data_homogeneous(iom, blockid=0, timerange=None):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    BF = BlockFactory()

    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)
    if timerange is not None:
        if len(timerange) == 1:
            I = (timegrid == timerange)
        else:
            I = ((timegrid >= timerange[0]) & (timegrid <= timerange[1]))
        if any(I):
            timegrid = timegrid[I]
        else:
            raise ValueError("No valid timestep remains!")

    # Basis shapes
    bsdescr = iom.load_wavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in bsdescr.items():
        BS[ahash] = BF.create_basis_shape(descr)

    # Plot the coefficients for all timesteps
    for j, step in enumerate(timegrid):
        allhashes, allcoeffs = iom.load_wavepacket_coefficients(timestep=step, blockid=blockid, get_hashes=True)

        k = []
        ck = []
        for ahash, coeffs in zip(allhashes, allcoeffs):
            bs = BS[int(ahash)]
            ki = array([ bs[node] for node in bs.get_node_iterator(mode="mag")])
            ck.append(coeffs[ki])
            ki.sort()
            k.append(ki)

        dt = parameters["dt"] if "dt" in parameters else None
        plot_coefficients(k, ck, step, dt, blockid=blockid)


def read_data_inhomogeneous(iom, blockid=0, timerange=None):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    BF = BlockFactory()

    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
    if timerange is not None:
        if len(timerange) == 1:
            I = (timegrid == timerange)
        else:
            I = ((timegrid >= timerange[0]) & (timegrid <= timerange[1]))
        if any(I):
            timegrid = timegrid[I]
        else:
            raise ValueError("No valid timestep remains!")

    # Basis shapes
    bsdescr = iom.load_inhomogwavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in bsdescr.items():
        BS[ahash] = BF.create_basis_shape(descr)

    # Plot the coefficients for all timesteps
    for j, step in enumerate(timegrid):
        allhashes, allcoeffs = iom.load_inhomogwavepacket_coefficients(timestep=step, blockid=blockid, get_hashes=True)

        k = []
        ck = []
        for ahash, coeffs in zip(allhashes, allcoeffs):
            bs = BS[int(ahash)]
            ki = array([ bs[node] for node in bs.get_node_iterator(mode="mag")])
            ck.append(coeffs[ki])
            ki.sort()
            k.append(ki)

        dt = parameters["dt"] if "dt" in parameters else None
        plot_coefficients(k, ck, step, dt, blockid=blockid)


def plot_coefficients(k, c, step, dt, blockid=0, path='.'):
    """
    :param parameters: A :py:class:`ParameterProvider` instance.
    :param timegrid: The timegrid that belongs to the coefficient values.
    :param coeffs: The coefficient values.
    :param imgsize: The size of the plot. For a large number of plotted
                    coefficients, we might have to increase this value.
    """
    print("Plotting the coefficients of data block '%s' at timestep %d" % (blockid,step))

    N = len(k)

    fig = figure()

    for n in range(N):
        ax = fig.add_subplot(N,1,n+1)

        stemcf(k[n], angle(c[n]), abs(c[n]))

        # axis formatting:
        m = max(abs(c[n]))
        ax.set_xlim(-1, max(k[n])+1)
        ax.set_ylim(-0.1*m, 1.1*m)

        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$c_k$")

    if dt is not None:
        fig.suptitle(r"Coefficients $c_k$ at time $t=%f$" % (step*dt))
    else:
        fig.suptitle(r"Coefficients $c_k$")

    fig.savefig(os.path.join(path, "wavepacket_coefficients_block_%s_timestep_%07d.png" % (blockid, step)))
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

    parser.add_argument("-t", "--timerange",
                        type = int,
                        help = "Plot only timestep(s) in this range",
                        nargs = "+",
                        default = None)

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
    if not "all" in args.blockid:
        blockids = [ bid for bid in args.blockid if bid in blockids ]

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting wavepacket coefficients in data block '%s'" % blockid)

        if iom.has_wavepacket(blockid=blockid):
            read_data_homogeneous(iom, blockid=blockid, timerange=args.timerange, path=resultspath)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            read_data_inhomogeneous(iom, blockid=blockid, timerange=args.timerange, path=resultspath)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '%s'" % blockid)

    iom.finalize()
