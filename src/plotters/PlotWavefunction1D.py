"""The WaveBlocks Project

Plot the wavefunctions probability densities
for one-dimensional wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import angle, conj, real, imag, squeeze
from matplotlib.pyplot import figure, close

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import IOManager
from WaveBlocksND.Plot import plotcf
from WaveBlocksND import GlobalDefaults as GLD


def plot_frames(PP, iom, blockid=0, view=None, plotphase=True, plotcomponents=False, plotabssqr=False, load=True, gridblockid=None, imgsize=(12,9)):
    """Plot the wave function for a series of timesteps.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param view: The aspect ratio.
    :param plotphase: Whether to plot the complex phase. (slow)
    :param plotcomponents: Whether to plot the real/imaginary parts..
    :param plotabssqr: Whether to plot the absolute value squared.
    """
    parameters = iom.load_parameters()

    if not parameters["dimension"] == 1:
        print("No wavefunction of one space dimensions, silent return!")
        return

    if PP is None:
        PP = parameters

    if load is True:
        if gridblockid is None:
            gridblockid = blockid
        print("Loading grid data from datablock '%s'" % gridblockid)
        G = iom.load_grid(blockid=gridblockid)
        G = G.reshape((1, -1))
        grid = real(squeeze(G))
    else:
        print("Creating new grid")
        G = BlockFactory().create_grid(PP)
        grid = real(squeeze(G.get_nodes(flat=True)))

    # View
    if view[0] is None:
        view[0] = grid.min()
    if view[1] is None:
        view[1] = grid.max()

    timegrid = iom.load_wavefunction_timegrid(blockid=blockid)

    for step in timegrid:
        print(" Plotting frame of timestep # %d" % step)

        wave = iom.load_wavefunction(blockid=blockid, timestep=step)
        values = [ wave[j,...] for j in xrange(parameters["ncomponents"]) ]

        # Plot the probability densities projected to the eigenbasis
        fig = figure(figsize=imgsize)

        for index, component in enumerate(values):
            ax = fig.add_subplot(parameters["ncomponents"],1,index+1)
            ax.ticklabel_format(style="sci", scilimits=(0,0), axis="y")

            if plotcomponents is True:
                ax.plot(grid, real(component))
                ax.plot(grid, imag(component))
                ax.set_ylabel(r"$\Re \varphi_"+str(index)+r", \Im \varphi_"+str(index)+r"$")
            if plotabssqr is True:
                ax.plot(grid, real(component*conj(component)))
                ax.set_ylabel(r"$\langle \varphi_"+str(index)+r"| \varphi_"+str(index)+r"\rangle$")
            if plotphase is True:
                plotcf(grid, angle(component), real(component*conj(component)))
                ax.set_ylabel(r"$\langle \varphi_"+str(index)+r"| \varphi_"+str(index)+r"\rangle$")

            ax.set_xlabel(r"$x$")

            # Set the aspect window
            ax.set_xlim(view[:2])
            ax.set_ylim(view[2:])

        if parameters.has_key("dt"):
            fig.suptitle(r"$\Psi$ at time $"+str(step*parameters["dt"])+r"$")
        else:
            fig.suptitle(r"$\Psi$")

        fig.savefig("wavefunction_block_%s_timestep_%07d.png" % (blockid, step))
        close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-p", "--paramfile",
                        type = str,
                        help = "The configuration parameter file",
                        nargs = "?",
                        default = None)

    parser.add_argument("-b", "--blockid",
                        help = "The data block to handle",
                        nargs = "*",
                        default = [0])

    parser.add_argument("-x", "--xrange",
                        type = float,
                        help = "The plot range on the x-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-y", "--yrange",
                        type = float,
                        help = "The plot range on the y-axis",
                        nargs = 2,
                        default = [0, 3.5])

    parser.add_argument("--plotphase",
                        action = "store_false",
                        help = "Plot the complex phase (slow)")

    parser.add_argument("--plotcomponents",
                        action = "store_true",
                        help = "Plot the real/imaginary parts.")

    parser.add_argument("--plotabssqr",
                        action = "store_true",
                        help = "Plot the absolute value squared.")

    args = parser.parse_args()

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=args.datafile)

    # Read file with parameter data for grid
    if args.paramfile:
        PL = ParameterLoader()
        PP = PL.load_from_file(args.paramfile)
    else:
        PP = None

    # The axes rectangle that is plotted
    view = args.xrange + args.yrange

    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        print("Plotting frames of data block '%s'" % blockid)
        # See if we have wavefunction values
        if iom.has_wavefunction(blockid=blockid):
            plot_frames(PP, iom,
                        blockid=blockid,
                        view=view,
                        plotphase=args.plotphase,
                        plotcomponents=args.plotcomponents,
                        plotabssqr=args.plotabssqr,
                        load=False)
        else:
            print("Warning: Not plotting any wavefunctions in block '%s'!" % blockid)

    iom.finalize()
