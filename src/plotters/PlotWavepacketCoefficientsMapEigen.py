"""The WaveBlocks Project

Plot the evolution of the coefficients $c_i$ of each component
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

Note that this script will fail to produce correct results if the
basis shapes are adaptive and their mappings mu incompatible!

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import abs, angle, transpose, hstack
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory
from WaveBlocksND import BasisTransformationHAWP
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Plot import plotcm
import GraphicsDefaults as GD


def read_all_datablocks(iom):
    r"""Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    parameters = iom.load_parameters()

    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_wavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_homogeneous(iom, blockid=blockid), blockid=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_coefficients(parameters, read_data_inhomogeneous(iom, blockid=blockid), blockid=blockid)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '%s'" % blockid)


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)
    dt = parameters["dt"] if parameters.has_key("dt") else 1.0
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
    CI = [ [] for i in xrange(HAWP.get_number_components()) ]

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
    dt = parameters["dt"] if parameters.has_key("dt") else 1.0
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
    CI = [ [] for i in xrange(HAWP.get_number_components()) ]

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
        fig.savefig("wavepacket_coefficients_map_eigen_abs_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        ax.matshow(angle(coeff))
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_eigen_angle_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=False, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_eigen_cm_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
        close(fig)

        fig = figure(figsize=imgsize)
        ax = fig.gca()
        plotcm(coeff, angle(coeff), abs(coeff), darken=True, axes=ax)
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$t$")
        ax.set_title(r"Coefficients $c_k^{%d}$" % jndex)
        fig.savefig("wavepacket_coefficients_map_eigen_cm_darken_block_"+str(blockid)+"_component_"+str(jndex)+GD.output_format)
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
                        default = ['0'])

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=args.datafile)

    # Read the data and plot it, one plot for each data block.
    read_all_datablocks(iom)

    iom.finalize()
