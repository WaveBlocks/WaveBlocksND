#!/usr/bin/env python
"""The WaveBlocks Project

Plot the wavepackets probability densities
for two-dimensional wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
from numpy import real
from matplotlib.pyplot import figure, close

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import BasisTransformationHAWP
from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Plot import plotcf2d


def plot_frames(PP, iom, blockid=0, load=False, eigentransform=False, timerange=None, view=None, path='.'):
    """Plot the wavepacket for a series of timesteps.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    """
    parameters = iom.load_parameters()
    BF = BlockFactory()

    if not parameters["dimension"] == 2:
        print("No two-dimensional wavepacket, silent return!")
        return

    if PP is None:
        PP = parameters

    if load is True:
        # TODO: Implement reshaping
        raise NotImplementedError("Loading of 2D grids is not implemented")
    else:
        G = BF.create_grid(PP)

    if eigentransform:
        V = BF.create_potential(parameters)
        BT = BasisTransformationHAWP(V)

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

    u, v = G.get_axes()
    u = real(u.reshape(-1))
    v = real(v.reshape(-1))

    # View
    if view is not None:
        if view[0] is None:
            view[0] = u.min()
        if view[1] is None:
            view[1] = u.max()
        if view[2] is None:
            view[2] = v.min()
        if view[3] is None:
            view[3] = v.max()

    for step in timegrid:
        print(" Plotting frame of timestep # {}".format(step))

        HAWP = iom.load_wavepacket(step, blockid=blockid)
        N = HAWP.get_number_components()

        # Transform the values to the eigenbasis
        if eigentransform:
            BT.transform_to_eigen(HAWP)

        psi = HAWP.evaluate_at(G.get_nodes(), prefactor=True, component=0)

        # Plot
        fig = figure()

        for level in range(N):
            z = psi[level]
            z = z.reshape(G.get_number_nodes())

            fig.add_subplot(N, 1, level + 1)
            plotcf2d(u, v, z, darken=0.3, limits=view)

        fig.savefig(os.path.join(path, "wavepacket_block_%s_level_%d_timestep_%07d.png" % (blockid, level, step)))
        close(fig)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datafile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = GLD.file_resultdatafile)

    parser.add_argument("-p", "--parametersfile",
                        type = str,
                        help = "The simulation data file",
                        nargs = "?",
                        default = None)

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

    parser.add_argument("-et", "--eigentransform",
                        action = "store_true",
                        help = "Transform the data into the eigenbasis before plotting")

    parser.add_argument("-x", "--xrange",
                        type = float,
                        help = "The plot range on the x-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-y", "--yrange",
                        type = float,
                        help = "The plot range on the y-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-t", "--timerange",
                        type = int,
                        help = "Plot only timestep(s) in this range",
                        nargs = "+",
                        default = None)

    args = parser.parse_args()


    # File with the simulation data
    resultspath = os.path.abspath(args.resultspath)

    if not os.path.exists(resultspath):
        raise IOError("The results path does not exist: {}".format(args.resultspath))

    datafile = os.path.abspath(os.path.join(args.resultspath, args.datafile))
    parametersfile = os.path.abspath(os.path.join(args.resultspath, args.parametersfile))

    # Read file with simulation data
    iom = IOManager()
    iom.open_file(filename=datafile)

    # Read file with parameter data for grid
    if args.parametersfile:
        PL = ParameterLoader()
        PP = PL.load_from_file(parametersfile)
    else:
        PP = None

    # Which blocks to handle
    blockids = iom.get_block_ids()
    if "all" not in args.blockid:
        blockids = [bid for bid in args.blockid if bid in blockids]

    # The axes rectangle that is plotted
    view = args.xrange + args.yrange

    # Iterate over all blocks
    for blockid in blockids:
        print("Plotting frames of data block '{}'".format(blockid))
        # See if we have wavepacket values
        if iom.has_wavepacket(blockid=blockid):
            plot_frames(PP, iom, blockid=blockid,
                        eigentransform=args.eigentransform,
                        timerange=args.timerange,
                        view=view,
                        path=resultspath)
        else:
            print("Warning: Not plotting any wavepackets in block '{}'".format(blockid))

    iom.finalize()
