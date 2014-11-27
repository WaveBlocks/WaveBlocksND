"""The WaveBlocks Project

Plot the coefficients $c_i$ of each component of a homogeneous or
inhomogeneous Hagedorn wavepacket for all timesteps during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import abs, angle, array
from matplotlib.pyplot import figure, close

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory
from WaveBlocksND.Plot import stemcf
from WaveBlocksND import GlobalDefaults as GLD
import GraphicsDefaults as GD


def read_all_datablocks(iom):
    r"""Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_wavepacket(blockid=blockid):
            read_data_homogeneous(iom, blockid=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            read_data_inhomogeneous(iom, blockid=blockid)
        else:
            print("Warning: Not plotting wavepacket coefficients in block '"+str(blockid)+"'!")


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    BF = BlockFactory()

    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)

    # Basis shapes
    bsdescr = iom.load_wavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in bsdescr.iteritems():
        BS[ahash] = BF.create_basis_shape(descr)

    # Plot the coefficients for all timesteps
    for j, step in enumerate(timegrid):
        hashes, coeffs = iom.load_wavepacket_coefficients(timestep=step, blockid=blockid, get_hashes=True)

        k = []

        for i in xrange(parameters["ncomponents"]):
            bs = BS[int(hashes[i])]
            ki = array([bs[node] for node in bs.get_node_iterator()])
            k.append(ki)

        dt = parameters["dt"] if parameters.has_key("dt") else None
        plot_coefficients(k, coeffs, step, dt, index=blockid)


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    BF = BlockFactory()

    parameters = iom.load_parameters()
    timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)

    # Basis shapes
    bsdescr = iom.load_inhomogwavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in bsdescr.iteritems():
        BS[ahash] = BF.create_basis_shape(descr)

    # Plot the coefficients for all timesteps
    for j, step in enumerate(timegrid):
        hashes, coeffs = iom.load_inhomogwavepacket_coefficients(timestep=step, blockid=blockid, get_hashes=True)

        k = []

        for i in xrange(parameters["ncomponents"]):
            bs = BS[int(hashes[i])]
            ki = array([bs[node] for node in bs.get_node_iterator()])
            k.append(ki)

        dt = parameters["dt"] if parameters.has_key("dt") else None
        plot_coefficients(k, coeffs, step, dt, index=blockid)


def plot_coefficients(k, c, step, dt, index=0):
    """
    :param parameters: A :py:class:`ParameterProvider` instance.
    :param timegrid: The timegrid that belongs to the coefficient values.
    :param coeffs: The coefficient values.
    :param imgsize: The size of the plot. For a large number of plotted
                    coefficients, we might have to increase this value.
    """
    print("Plotting the coefficients of data block '"+str(index)+"' at timestep "+str(step))

    N = len(k)

    fig = figure()

    for n in xrange(N):
        ax = fig.add_subplot(N,1,n+1)

        stemcf(k[n], angle(c[n]), abs(c[n]))

        # axis formatting:
        m = max(abs(c[n]))
        ax.set_xlim(-1, max(k[n])+1)
        ax.set_ylim(-0.1*m, 1.1*m)

        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$c_k$")

    if dt is not None:
        fig.suptitle(r"Coefficients $c_k$ at time $t="+str(step*dt)+r"$")
    else:
        fig.suptitle(r"Coefficients $c_k$")

    fig.savefig("wavepacket_coefficients_block"+str(index)+"_timestep_"+(5-len(str(step)))*"0"+str(step)+GD.output_format)
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
