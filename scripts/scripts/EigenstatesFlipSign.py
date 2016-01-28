#!/usr/bin/env python
"""The WaveBlocks Project

Flip the signs of the coefficients of some eigenstate wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse

from WaveBlocksND import IOManager


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--datafile",
                    type = str,
                    help = "The data file",
                    nargs = "?",
                    default = "eigenstates.hdf5")

parser.add_argument("-b", "--blockid",
                    help = "The data block to handle",
                    nargs = "*",
                    default = [0])

args = parser.parse_args()

# Read file
iom = IOManager()
iom.open_file(filename=args.datafile)

# Which blocks to handle
blockids = iom.get_block_ids()
if not "all" in args.blockid:
    blockids = [ bid for bid in args.blockid if bid in blockids ]

# Iterate over all blocks
for blockid in blockids:
    if iom.has_wavepacket(blockid=blockid):
        # Ugly hack using raw hdf5 data access
        path = "/datablock_"+str(blockid)+"/wavepacket/coefficients/c_0"
        iom._srf[path][:] *= -1.0

iom.finalize()
