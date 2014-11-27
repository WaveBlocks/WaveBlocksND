"""The WaveBlocks Project

Plot the wavefunctions probability densities
for two-dimensional wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import angle, real
from mayavi import mlab

from WaveBlocksND import ParameterLoader
from WaveBlocksND import BlockFactory
from WaveBlocksND import WaveFunction
from WaveBlocksND import BasisTransformationWF
from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GLD
from WaveBlocksND.Plot3D import surfcf


def plot_frames(PP, iom, blockid=0, load=False, tte=False, view=None, interactive=False, sparsify=10):
    """Plot the wave function for a series of timesteps.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param view: The aspect ratio.
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

    if tte:
        V = BlockFactory().create_potential(parameters)
        BT = BasisTransformationWF(V)
        BT.set_grid(G)

    WF = WaveFunction(parameters)
    WF.set_grid(G)
    N = WF.get_number_components()

    timegrid = iom.load_wavefunction_timegrid(blockid=blockid)

    u, v = G.get_nodes(split=True, flat=False)
    u = real(u[::sparsify,::sparsify])
    v = real(v[::sparsify,::sparsify])

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
        print(" Plotting frame of timestep # %d" % step)

        # Load the data
        wave = iom.load_wavefunction(blockid=blockid, timestep=step)
        values = [ wave[j,...] for j in xrange(parameters["ncomponents"]) ]
        WF.set_values(values)

        # Transform the values to the eigenbasis
        if tte:
            BT.transform_to_eigen(WF)

        Psi = WF.get_values()

        for level in xrange(N):
            # Wavefunction data
            z = Psi[level]
            z = z.reshape(G.get_number_nodes())[::sparsify,::sparsify]

            # View
            if view is not None:
                if view[4] is None:
                    view[4] = 0.0
                if view[5] is None:
                    view[5] = 1.1*abs(z).max()

            # Plot
            if not interactive:
                mlab.options.offscreen = True

            fig = mlab.figure(size=(800,700))

            surfcf(u, v, angle(z), abs(z), view=view)

            mlab.draw()
            if interactive:
                mlab.show()
            else:
                mlab.savefig("wavefunction_block_%s_level_%d_timestep_%07d.png" % (blockid, level, step))
                mlab.close(fig)




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
                        default = ['0'])

    parser.add_argument("-tte", "--transformtoeigen",
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

    parser.add_argument("-z", "--zrange",
                        type = float,
                        help = "The plot range on the z-axis",
                        nargs = 2,
                        default = [None, None])

    parser.add_argument("-s", "--sparsify",
                        type = int,
                        help = "Plot only every s-th point",
                        default = 10)

    parser.add_argument("-i", "--interactive",
                        action = "store_true",
                        help = "Hold the plot open for interactive manipulation")

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
    view = args.xrange + args.yrange + args.zrange

    # Iterate over all blocks and plot their data
    for blockid in iom.get_block_ids():
        print("Plotting frames of data block '%s'" % blockid)
        # See if we have wavefunction values
        if iom.has_wavefunction(blockid=blockid):
            plot_frames(PP, iom, blockid=blockid,
                        tte=args.transformtoeigen,
                        sparsify=args.sparsify,
                        view=view,
                        interactive=args.interactive)
        else:
            print("Warning: Not plotting any wavefunctions in block '%s'" % blockid)

    iom.finalize()
