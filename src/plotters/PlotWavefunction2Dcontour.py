"""The WaveBlocks Project

Plot the wavefunctions probability densities
for two-dimensional wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from numpy import real
from matplotlib.pyplot import *

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import WaveFunction
from WaveBlocksND import BasisTransformationWF
from WaveBlocksND import IOManager

from WaveBlocksND.Plot import plotcm


def plot_frames(PP, iom, blockid=0, load=False):
    r"""
    """
    parameters = iom.load_parameters()

    if not parameters["dimension"] == 2:
        print("No wavefunction of two space dimensions, silent return!")
        return

    if PP is None:
        PP = parameters

    if load is True:
        # TODO: Implement reshaping
        raise NotImplementedError("Loading of 2D grids is not implemented")
        #G = iom.load_grid(blockid=blockid)
        #G = grid.reshape((1, -1))
    else:
        G = BlockFactory().create_grid(PP)

    V = BlockFactory().create_potential(parameters)

    WF = WaveFunction(parameters)
    WF.set_grid(G)

    BT = BasisTransformationWF(V)
    BT.set_grid(G)

    timegrid = iom.load_wavefunction_timegrid(blockid=blockid)

    u, v = G.get_nodes(split=True, flat=False)
    u = real(u)
    v = real(v)

    N = WF.get_number_components()

    for step in timegrid:
        print(" Plotting frame of timestep # " + str(step))

        wave = iom.load_wavefunction(blockid=blockid, timestep=step)
        values = [ wave[j,...] for j in xrange(parameters["ncomponents"]) ]

        WF.set_values(values)

        # Transform the values to the eigenbasis
        # TODO: improve this:
        if parameters["algorithm"] == "fourier":
            BT.transform_to_eigen(WF)
        else:
            pass

        Psi = WF.get_values()

        fig = figure()

        for level in xrange(N):
            z = Psi[level]
            z = z.reshape(G.get_number_nodes())

            subplot(N,1,level+1)
            plotcm(z, darken=0.3)

        savefig("wavefunction_level_"+str(level)+"_timestep_"+(5-len(str(step)))*"0"+str(step)+".png")
        close(fig)

    print(" Plotting frames finished")



if __name__ == "__main__":
    iom = IOManager()
    PL = ParameterLoader()

    # Read file with simulation data
    try:
        iom.open_file(filename=sys.argv[1])
    except IndexError:
        iom.open_file()

    # Read file with parameter data for grid
    try:
        PP = PL.load_from_file(sys.argv[2])
    except IndexError:
        PP = None

    # The axes rectangle that is plotted
    #view = [-3.5, 3.5, -0.1, 3.5]

    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        print("Plotting frames of data block '"+str(blockid)+"'")
        # See if we have wavefunction values
        if iom.has_wavefunction(blockid=blockid):
            plot_frames(PP, iom, blockid=blockid)
        else:
            print("Warning: Not plotting any wavefunctions in block '"+str(blockid)+"'!")

    iom.finalize()
