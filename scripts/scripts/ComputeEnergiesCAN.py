#!/usr/bin/env python
"""The WaveBlocks Project

Compute the energies of the different wavepackets or wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os

from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GD

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--datafile",
                    type = str,
                    help = "The simulation data file.",
                    nargs = "?",
                    default = GD.file_resultdatafile)

parser.add_argument("-b", "--blockid",
                    type = str,
                    help = "The data block to handle.",
                    nargs = "*",
                    default = ["all"])

parser.add_argument("-r", "--resultspath",
                    type = str,
                    help = "Path where to put the results.",
                    nargs = "?",
                    default = '.')

parser.add_argument("-et", "--eigentransform",
                    help = "Transform the data into the eigenbasis before computing norms.",
                    action = "store_true")

# TODO: Filter type of objects
# parser.add_argument("-t", "--type",
#                     help = "The type of objects to consider.",
#                     type = str,
#                     default = "all")

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

# Iterate over all blocks
for blockid in blockids:
    print("Computing the energies in data block '%s'" % blockid)

    if iom.has_energy(blockid=blockid):
        print("Datablock '%s' already contains energy data, silent skip." % blockid)
        continue

    # NOTE: Add new algorithms here

    if iom.has_wavepacket(blockid=blockid):
        from WaveBlocksND.Interface import EnergiesWavepacket
        EnergiesWavepacket.compute_energy_hawp(iom, blockid=blockid, eigentrafo=args.eigentransform, iseigen=False)
    elif iom.has_wavefunction(blockid=blockid):
        from WaveBlocksND.Interface import EnergiesWavefunction
        EnergiesWavefunction.compute_energy(iom, blockid=blockid, eigentrafo=args.eigentransform, iseigen=False)
    elif iom.has_inhomogwavepacket(blockid=blockid):
        from WaveBlocksND.Interface import EnergiesWavepacket
        EnergiesWavepacket.compute_energy_inhawp(iom, blockid=blockid, eigentrafo=args.eigentransform, iseigen=False)
    else:
        print("Warning: Not computing any energies in block '%s'!" % blockid)

iom.finalize()
