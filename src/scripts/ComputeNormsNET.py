"""The WaveBlocks Project

Compute the norms of the different wavepackets or wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014 R. Bourquin
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

parser.add_argument("-et", "--eigentransform",
                    help = "Transform the data into the eigenbasis before computing norms",
                    action = "store_false")

# TODO: Filter type of objects
# parser.add_argument("-t", "--type",
#                     help = "The type of objects to consider",
#                     type = str,
#                     default = "all")

args = parser.parse_args()

# Read file with simulation data
iom = IOManager()
iom.open_file(filename=args.datafile)

# Which blocks to handle
blockids = iom.get_block_ids()
if not "all" in args.blockid:
    blockids = [ bid for bid in args.blockid if bid in blockids ]

# Iterate over all blocks
for blockid in blockids:
    print("Computing the norms in data block '%s'" % blockid)

    if iom.has_norm(blockid=blockid):
        print("Datablock '%s' already contains norm data, silent skip" % blockid)
        continue

    # NOTE: Add new algorithms here

    if iom.has_wavepacket(blockid=blockid):
        import NormWavepacket
        NormWavepacket.compute_norm_hawp(iom, blockid=blockid, eigentrafo=args.eigentransform)
    elif iom.has_wavefunction(blockid=blockid):
        import NormWavefunction
        NormWavefunction.compute_norm(iom, blockid=blockid, eigentrafo=args.eigentransform)
    elif iom.has_inhomogwavepacket(blockid=blockid):
        import NormWavepacket
        NormWavepacket.compute_norm_inhawp(iom, blockid=blockid, eigentrafo=args.eigentransform)
    else:
        print("Warning: Not computing any norm in block '%s'!" % blockid)

iom.finalize()
