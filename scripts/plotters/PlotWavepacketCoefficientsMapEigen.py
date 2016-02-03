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
from numpy import abs, angle, transpose, hstack
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory
from WaveBlocksND import BasisTransformationHAWP
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Plot import plotcm
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

    # The potential used
    Potential = BlockFactory().create_potential(parameters)

    # Basis transformator
    BT = BasisTransformationHAWP(Potential)

    # Initialize a Hagedorn wavepacket with the data
    descr = iom.load_wavepacket_description(blockid=blockid)
    HAWP = BlockFactory().create_wavepacket(descr)

    BT.set_matrix_builder(HAWP.get_innerproduct())

    # Store the resulting coefficients here
    CI = [ [] for i in range(HAWP.get_number_components()) ]

    # Iterate over all timesteps, this is an *expensive* transformation
    for i, step in enumerate(timegrid):
        print(" Computing eigentransformation at timestep %d" % step)
        # Retrieve simulation data
        HAWP = iom.load_wavepacket(timestep=step, blockid=blockid)

        # Transform to the eigenbasis.
        BT.transform_to_eigen(HAWP)

        for index, item in enumerate(HAWP.get_coefficients()):
            CI[index].append(item)

    CI = [ transpose(hstack(item)) for item in CI ]

    return time, CI


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if "dt" in parameters else 1.0
    time = timegrid * dt

    # The potential used
    Potential = BlockFactory().create_potential(parameters)

    # Basis transformator
    BT = BasisTransformationHAWP(Potential)

    # Initialize a Hagedorn wavepacket with the data
    descr = iom.load_inhomogwavepacket_description(blockid=blockid)
    HAWP = BlockFactory().create_wavepacket(descr)

    BT.set_matrix_builder(HAWP.get_quadrature())

    # Store the resulting coefficients here
    CI = [ [] for i in range(HAWP.get_number_components()) ]

    # Iterate over all timesteps, this is an *expensive* transformation
    for i, step in enumerate(timegrid):
        print(" Computing eigentransformation at timestep %d" % step)
        # Retrieve simulation data
        HAWP = iom.load_inhomogwavepacket(timestep=step, blockid=blockid)

        # Transform to the eigenbasis.
        BT.transform_to_eigen(HAWP)

        for index, item in enumerate(HAWP.get_coefficients()):
            CI[index].append(item)

    CI = [ transpose(hstack(item)) for item in CI ]

    return time, CI


def plot_coefficients(parameters, data, blockid=0, imgsize=(10,20), path='.'):
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
        fig.savefig(os.path.join(path, "wavepacket_coefficients_map_eigen_abs_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format))
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        ax.matshow(angle(coeff))
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig(os.path.join(path, "wavepacket_coefficients_map_eigen_angle_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format))
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=False, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig(os.path.join(path, "wavepacket_coefficients_map_eigen_cm_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format))
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=True, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig(os.path.join(path, "wavepacket_coefficients_map_eigen_cm_darken_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format))
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

    parameters = iom.load_parameters()

    # Which blocks to handle
    blockids = iom.get_block_ids()
    if not "all" in args.blockid:
        blockids = [ bid for bid in args.blockid if bid in blockids ]

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting wavepacket coefficients in data block '%s'" % blockid)

        # NOTE: Add new algorithms here

        if iom.has_wavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_homogeneous(iom, blockid=blockid), blockid=blockid, path=resultspath)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_inhomogeneous(iom, blockid=blockid), blockid=blockid, path=resultspath)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '%s'" % blockid)

    iom.finalize()
