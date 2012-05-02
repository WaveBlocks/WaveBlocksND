r"""The WaveBlocks Project

Plot the evolution of the parameters :math:`Pi_i = (q, p, Q, P, S)`
of a homogeneous or inhomogeneous Hagedorn wavepacket during the
time propagation.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from numpy import real, imag, abs, squeeze
from numpy.linalg import norm, det
from matplotlib.mlab import amap
from matplotlib.pyplot import *

from WaveBlocksND import ComplexMath
from WaveBlocksND import IOManager

import GraphicsDefaults as GD


def read_all_datablocks(iom):
    r"""Read the data from all blocks that contain any usable data.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        if iom.has_wavepacket(blockid=blockid):
            plot_parameters(read_data_homogeneous(iom, blockid=blockid), index=blockid)
        elif iom.has_inhomogwavepacket(blockid=blockid):
            plot_parameters(read_data_inhomogeneous(iom, blockid=blockid), index=blockid)
        else:
            print("Warning: Not plotting wavepacket parameters in block '"+str(blockid)+"'!")


def read_data_homogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
    timegrid = iom.load_wavepacket_timegrid(blockid=blockid)
    time = timegrid * parameters["dt"]

    Pi = iom.load_wavepacket_parameters(blockid=blockid)

    # The Dimension D, we know that q has shapw (#timesteps, D, 1)
    D = Pi[0].shape[1]
    if not D == 2:
        print("Warning: Trajectory plotting implemented only for 2D wavepackets")
        return

    Pi = map(squeeze, Pi)
    qhist, phist, Qhist, Phist, Shist = Pi
    return (time, [qhist], [phist], [Qhist], [Phist], [Shist])


def read_data_inhomogeneous(iom, blockid=0):
    r"""
    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param blockid: The data block from which the values are read.
    """
    parameters = iom.load_parameters()
#     timegrid = iom.load_inhomogwavepacket_timegrid(blockid=blockid)
#     time = timegrid * parameters["dt"]

#     Pi = iom.load_inhomogwavepacket_parameters(blockid=blockid)

#     # Number of components
#     N = parameters["ncomponents"]

#     Phist = [ Pi[i][:,0] for i in xrange(N) ]
#     Qhist = [ Pi[i][:,1] for i in xrange(N) ]
#     Shist = [ Pi[i][:,2] for i in xrange(N) ]
#     phist = [ Pi[i][:,3] for i in xrange(N) ]
#     qhist = [ Pi[i][:,4] for i in xrange(N) ]

#    return (time, Phist, Qhist, Shist, phist, qhist)


def plot_parameters(data, index=0):
    r"""Plot the data parameters :math:`(q,p,Q,P,S)` over time.
    For each new `index` we start a new figure. This allows plotting
    several time evolutions to the same figure
    """
    print("Plotting the parameters of data block '"+str(index)+"'")

    timegrid, qhist, phist, Qhist, Phist, Shist = data

    # Plot the 2D trajectory of the parameters q and p
    fig = figure()
    ax = fig.gca()
    for item in qhist:
        ax.plot(item[:,0], item[:,1], "-o", label=r"Trajectory of $q$")
    ax.grid(True)
    ax.set_title(r"Trajectory of $q$")
    fig.savefig("wavepacket_parameters_trajectoryq_block"+str(index)+GD.output_format)
    close(fig)

    # Plot the 2D trajectory of the parameters q and p
    fig = figure()
    ax = fig.gca()
    for item in phist:
        ax.plot(item[:,0], item[:,1], "-o", label=r"Trajectory of $p$")
    ax.grid(True)
    ax.set_title(r"Trajectory of $p$")
    fig.savefig("wavepacket_parameters_trajectoryp_block"+str(index)+GD.output_format)
    close(fig)




if __name__ == "__main__":
    iom = IOManager()

    # Read file with simulation data
    try:
        iom.open_file(filename=sys.argv[1])
    except IndexError:
        iom.open_file()

    # Read the data and plot it, one plot for each data block.
    read_all_datablocks(iom)

    iom.finalize()
