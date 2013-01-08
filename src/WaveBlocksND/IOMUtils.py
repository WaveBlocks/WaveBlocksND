"""The WaveBlocks Project

Various small utility functions.

@author: R. Bourquin
@copyright: Copyright (C) 2013 R. Bourquin
@license: Modified BSD License
"""

from BlockFactory import BlockFactory


def load_wavepacket(iom, timestep, blockid=1):

    if not iom.has_wavepacket(blockid=blockid):
        print("There is no wavepacket to load in block "+str(blockid))
        return

    BF = BlockFactory()

    wpd = iom.load_wavepacket_description(blockid=blockid)
    HAWP = BF.create_wavepacket(wpd)

    # Basis shapes
    BS_descr = iom.load_wavepacket_basisshapes(blockid=blockid)
    BS = {}
    for ahash, descr in BS_descr.iteritems():
        BS[ahash] = BF.create_basis_shape(descr)

    KEY = ("q","p","Q","P","S","adQ")
    # Retrieve simulation data
    params = iom.load_wavepacket_parameters(timestep=timestep, blockid=blockid, key=KEY)
    hashes, coeffs = iom.load_wavepacket_coefficients(timestep=timestep, get_hashes=True, blockid=blockid)
    # Configure the wavepacket
    HAWP.set_parameters(params, key=KEY)
    HAWP.set_basis_shapes([ BS[int(ha)] for ha in hashes ])
    HAWP.set_coefficients(coeffs)

    return HAWP
