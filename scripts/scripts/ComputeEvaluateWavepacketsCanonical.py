#!/usr/bin/env python
"""The WaveBlocks Project

Sample wavepackets at the nodes of a given grid and save
the results back to the given simulation data file.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse

from WaveBlocksND import IOManager
from WaveBlocksND import ParameterLoader
from WaveBlocksND import GlobalDefaults as GD

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--datafile",
                    type = str,
                    help = "The simulation data file",
                    nargs = "?",
                    default = GD.file_resultdatafile)

parser.add_argument("-b", "--blockid",
                    type = str,
                    help = "The data block to handle",
                    nargs = "*",
                    default = ["all"])

parser.add_argument("-p", "--paramfile",
                    type = str,
                    help = "The configuration parameter file",
                    nargs = "?",
                    default = None)

parser.add_argument("-et", "--eigentransform",
                    help = "Transform the data into the eigenbasis before computing norms",
                    action = "store_true")

args = parser.parse_args()

# Read file with simulation data
iom = IOManager()
iom.open_file(filename=args.datafile)

# Read the additional grid parameters
if args.paramfile:
    parametersfile = args.paramfile
    PA = ParameterLoader().load_from_file(parametersfile)
else:
    PA = None

# Which blocks to handle
blockids = iom.get_block_ids()
if not "all" in args.blockid:
    blockids = [ bid for bid in args.blockid if bid in blockids ]

# Iterate over all blocks
for blockid in blockids:
    print("Evaluating wavepackets in data block '%s'" % blockid)

    if iom.has_wavefunction(blockid=blockid):
        print("Datablock '%s' already contains wavefunction data, silent skip." % blockid)
        continue

    # NOTE: Add new algorithms here

    if iom.has_wavepacket(blockid=blockid):
        import EvaluateWavepackets
        EvaluateWavepackets.compute_evaluate_wavepackets(PA, iom, blockid=blockid, eigentrafo=args.eigentransform)
    elif iom.has_inhomogwavepacket(blockid=blockid):
        import EvaluateWavepacketsInhomog
        EvaluateWavepacketsInhomog.compute_evaluate_wavepackets(PA, iom, blockid=blockid, eigentrafo=args.eigentransform)
    else:
        print("Warning: Not evaluating any wavepackets in block '%s'!" % blockid)

iom.finalize()
