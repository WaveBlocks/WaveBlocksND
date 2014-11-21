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
if "all" in args.blockid:
    blocks_to_handle = iom.get_block_ids()
else:
    blocks_to_handle = map(int, args.blockid)

# Iterate over all blocks
for blockid in blocks_to_handle:
    try:
        # Ugly hack using raw hdf5 data access
        path = "/datablock_"+str(blockid)+"/wavepacket/coefficients/c_0"
        iom._srf[path][:] *= -1.0
    except KeyError:
        pass

iom.finalize()
