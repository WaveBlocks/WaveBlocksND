"""The WaveBlocks Project

This script loads and transforms data from a datafile
to make them suitable as initial values.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND import IOManager
from WaveBlocksND import BlockFactory


__all__ = ["load_from_file"]


def load_from_file(filepath, blockid=0, timestep=0):

    IOM = IOManager()
    IOM.open_file(filepath)

    # Check if we have data
    tg = IOM.load_wavepacket_timegrid(blockid=blockid)
    if not timestep in tg:
        raise ValueError("No data for timestep "+str(timestep))

    # Load data and assemble packet
    BF = BlockFactory()

    # Basis shapes
    BS_descr = IOM.load_wavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in BS_descr.iteritems():
        BS[ahash] = BF.create_basis_shape(descr)

    # Create a packet
    wpd = IOM.load_wavepacket_description(blockid=blockid)
    HAWP = BF.create_wavepacket(wpd)

    # Data
    ha, ci = IOM.load_wavepacket_coefficients(blockid=blockid, timestep=timestep, get_hashes=True)
    Pi = IOM.load_wavepacket_parameters(blockid=blockid, timestep=timestep)

    HAWP.set_parameters(Pi)
    HAWP.set_basis_shapes([ BS[int(h)] for h in ha ])
    HAWP.set_coefficients(ci)

    # Reformat data
    C = []

    for n in xrange(HAWP.get_number_components()):
        B = HAWP.get_basis_shapes(component=n)
        cn = HAWP.get_coefficients(component=n)
        l = []
        for i in xrange(B.get_basis_size()):
            l.append( (B[i], cn[i,0]) )
        C.append(l)

    return Pi, C
