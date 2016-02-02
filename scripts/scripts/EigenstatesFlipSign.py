#!/usr/bin/env python
"""The WaveBlocks Project

Flip the signs of the coefficients of some eigenstate wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os

from WaveBlocksND import IOManager

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--datafile",
                    type = str,
                    help = "The data file.",
                    nargs = "?",
                    default = "eigenstates.hdf5")

parser.add_argument("-b", "--blockid",
                    help = "The data block to handle.",
                    nargs = "*",
                    default = [0])

parser.add_argument("-r", "--resultspath",
                    type = str,
                    help = "Path where to put the results.",
                    nargs = "?",
                    default = '.')

args = parser.parse_args()


# File with the simulation data
resultspath = os.path.abspath(args.resultspath)

if not os.path.exists(resultspath):
    raise IOError("The results path does not exist: " + args.resultspath)

datafile = os.path.abspath(os.path.join(args.resultspath, args.datafile))

# Read file
iom = IOManager()
iom.open_file(filename=datafile)

# Which blocks to handle
blockids = iom.get_block_ids()
if not "all" in args.blockid:
    blockids = [ bid for bid in args.blockid if bid in blockids ]

# Iterate over all blocks
for blockid in blockids:
    if iom.has_wavepacket(blockid=blockid):
        # TODO: Ugly hack using raw hdf5 data access
        path = "/datablock_"+str(blockid)+"/wavepacket/coefficients/c_0"
        iom._srf[path][:] *= -1.0

iom.finalize()
