#!/usr/bin/env python
"""The WaveBlocks Project

Plot the evolution of the coefficients $c_i$ of each component
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

Note that this script will fail to produce correct results if the
basis shapes are adaptive and their mappings mu incompatible!

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real, imag, abs, angle
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory
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
    time = timegrid * dt

    hashes, coeffs = iom.load_wavepacket_coefficients(blockid=blockid, get_hashes=True)

    return time, coeffs, [h.reshape(-1) for h in hashes]


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if "dt" in parameters else 1.0
    time = timegrid * dt

    hashes, coeffs = iom.load_inhomogwavepacket_coefficients(blockid=blockid, get_hashes=True)

    return time, coeffs, [h.reshape(-1) for h in hashes]


def plot_coefficients(parameters, data, absang=False, index=0, reim=False, imgsize=(10,20), path='.'):
    """
    :param parameters: A :py:class:`ParameterProvider` instance.
    :param timegrid: The timegrid that belongs to the coefficient values.
    :param coeffs: The coefficient values.
    :param imgsize: The size of the plot. For a large number of plotted
                    coefficients, we might have to increase this value.
    """
    print("Plotting the coefficients of data block '%s'" % index)

    # Check if we have enough coefficients to plot
    timegrid, coeffs, hashes = data

    N = len(coeffs)

    # Reconstruct basis shapes
    BS = []
    for i in range(N):
        d = iom.load_wavepacket_basisshapes(the_hash=hashes[i][0])
        bs = BlockFactory().create_basis_shape(d)
        BS.append(bs)

    # Compute all multi indices available (union of all shapes)
    allvects = set([])

    for bs in BS:
        allvects.update(set([i for i in bs]))

    # Plot
    for vect in allvects:
        print(" Plotting coefficient "+str(vect)+" out of %d" % len(allvects))

        fig = figure()

        for level in range(N):
            j = BS[level][vect]
            if j is not None:
                ax = fig.add_subplot(N,1,level+1)

                if not reim:
                    ax.plot(timegrid, angle(coeffs[level][:,j]), label=r"$\arg c$")
                else:
                    ax.plot(timegrid, real(coeffs[level][:,j]), label=r"$\Re c$")
                    ax.plot(timegrid, imag(coeffs[level][:,j]), label=r"$\Im c$")

                ax.plot(timegrid, abs(coeffs[level][:,j]), "r", label=r"$|c|$")

                ax.grid(True)
                ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")
                ax.legend(loc='upper right')
                ax.set_xlabel(r"$t$")
                ax.set_ylabel(r"$c^{%d}$" % level)

        fig.suptitle(r"$\underline{k} = "+str(vect)+r"$")
        fig.savefig(os.path.join(path, "coefficient_k"+str(vect)+"_block"+str(index)+GD.output_format))
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

    parser.add_argument("--reim",
                        action = "store_true",
                        help = "Plot the real and imaginary parts")

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

    # Read the data and plot it, one plot for each data block.
    parameters = iom.load_parameters()

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting wavepacket coefficients in data block '%s'" % blockid)

        # NOTE: Add new algorithms here

        if iom.has_wavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_homogeneous(iom, blockid=blockid), index=blockid, reim=args.reim, path=resultspath)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_inhomogeneous(iom, blockid=blockid), index=blockid, reim=args.reim, path=resultspath)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '%s'" % blockid)

    iom.finalize()