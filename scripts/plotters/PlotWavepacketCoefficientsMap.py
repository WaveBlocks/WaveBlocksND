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
from numpy import abs, angle
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND.Plot import plotcm
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

    return time, coeffs


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

    return time, coeffs


def plot_coefficients(parameters, data, blockid=0, imgsize=(10,20)):
    """
    :param parameters: A :py:class:`ParameterProvider` instance.
    :param timegrid: The timegrid that belongs to the coefficient values.
    :param coeffs: The coefficient values.
    :param imgsize: The size of the plot. For a large number of plotted
                    coefficients, we might have to increase this value.
    """
    print("Plotting the coefficients of data block '%s'" % blockid)

    # Check if we have enough coefficients to plot
    timegrid, coeffs = data

    # Plot for each of the N levels
    for jndex, coeff in enumerate(coeffs):

        # Scale the image to roughly fit the data shape
        v, u = coeff.shape
        imags = (imgsize[0], int(10*v / (1.0*u)))
        if imags[1] < 10000:
            imgsize = imags

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        ax.matshow(abs(coeff))
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_abs_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        ax.matshow(angle(coeff))
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_angle_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=False, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_cm_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=True, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_cm_darken_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
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

    parameters = iom.load_parameters()

    # Which blocks to handle
    blockids = iom.get_block_ids()
    if not "all" in args.blockid:
        blockids = [ bid for bid in args.blockid if bid in blockids ]

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting wavepacket parameters in data block '%s'" % blockid)

        if iom.has_wavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_homogeneous(iom, blockid=blockid), blockid=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_inhomogeneous(iom, blockid=blockid), blockid=blockid)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '%s'" % blockid)

    iom.finalize()
