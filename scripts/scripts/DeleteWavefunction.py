#!/usr/bin/env python
"""The WaveBlocks Project

Delete wavefunction data from simulation result files.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse

from WaveBlocksND import IOManager
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

args = parser.parse_args()

# Read file with simulation data
iom = IOManager()
iom.open_file(filename=args.datafile)

# Which blocks to handle
blockids = iom.get_block_ids()
if "all" not in args.blockid:
    blockids = [bid for bid in args.blockid if bid in blockids]

# Iterate over all blocks
for blockid in blockids:
    if iom.has_wavefunction(blockid=blockid):
        print("Deleting grid and wavefunction data in block '%s'" % blockid)
        iom.delete_wavefunction(blockid=blockid)
        if iom.has_grid(blockid=blockid):
            iom.delete_grid(blockid=blockid)

iom.finalize()
