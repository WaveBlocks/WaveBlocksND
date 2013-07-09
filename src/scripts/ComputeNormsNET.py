"""The WaveBlocks Project

Compute the norms of the different wavepackets or wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013 R. Bourquin
@license: Modified BSD License
"""

import sys
import argparse

from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults as GD

parser = argparse.ArgumentParser()

parser.add_argument("simfile",
                    type = str,
                    help = "The simulation data file",
                    default = GD.file_resultdatafile)

parser.add_argument("-b", "--blockid",
                    help = "The data block to handle",
                    nargs = "*",
                    default = [0])

# TODO: Filter type of objects
# parser.add_argument("-t", "--type",
#                     help = "The type of objects to consider",
#                     type = str,
#                     default = "all")

args = parser.parse_args()

# Read file with simulation data
iom = IOManager()
iom.open_file(filename=args.simfile)

# Which blocks to handle
if "all" in args.blockid:
    blocks_to_handle = iom.get_block_ids()
else:
    blocks_to_handle = map(int, args.blockid)

# Iterate over all blocks
for blockid in blocks_to_handle:
    print("Computing the norms in data block '"+str(blockid)+"'")

    if iom.has_norm(blockid=blockid):
        print("Datablock '"+str(blockid)+"' already contains norm data, silent skip.")
        continue

    # TODO: Add new algorithms here

    # We test for a homogeneous wavepacket next
    if iom.has_wavepacket(blockid=blockid):
        import NormWavepacket
        NormWavepacket.compute_norm_hawp(iom, blockid=blockid, eigentrafo=False)
        # We have no wavepacket, then we try for a wavefunction
    elif iom.has_wavefunction(blockid=blockid):
        import NormWavefunction
        NormWavefunction.compute_norm(iom, blockid=blockid, eigentrafo=False)
    # We test for an inhomogeneous wavepacket next
    elif iom.has_inhomogwavepacket(blockid=blockid):
        import NormWavepacket
        NormWavepacket.compute_norm_inhawp(iom, blockid=blockid, eigentrafo=False)
    # If there is also no wavefunction, then there is nothing to compute the norm
    else:
        print("Warning: Not computing any norm in block '"+str(blockid)+"'!")

iom.finalize()
