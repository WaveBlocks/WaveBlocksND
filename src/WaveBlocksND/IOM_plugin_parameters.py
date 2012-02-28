"""The WaveBlocks Project

IOM plugin providing functions for handling simulation parameter data.

@author: R. Bourquin
@copyright: Copyright (C) 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import pickle

import ParameterProvider as ParameterProvider


def add_parameters(self, blockid="global"):
    """Add storage for the simulation parameters.
    """
    # Store the simulation parameters
    # We are only interested in the attributes of this data set
    # as they are used to store the simulation parameters.
    paset = self._srf[self._prefixb+str(blockid)].create_dataset("simulation_parameters", (1,1))


def delete_parameters(self, blockid="global"):
    """Remove the stored simulation parameters.
    """
    try:
        del self._srf[self._prefixb+str(blockid)+"/simulation_parameters"]
    except KeyError:
        pass


def has_parameters(self, blockid="global"):
    """Ask if the specified data block has the desired data tensor.
    """
    return "simulation_parameters" in self._srf[self._prefixb+str(blockid)].keys()


def save_parameters(self, parameters, blockid="global"):
    """Save the norm of wavefunctions or wavepackets.
    """
    paset = self._srf["/"+self._prefixb+str(blockid)+"/simulation_parameters"]

    for param, value in parameters:
        # Store all the values as pickled strings because hdf can
        # only store strings or ndarrays as attributes.
        paset.attrs[param] = pickle.dumps(value)


def load_parameters(self, blockid="global"):
    """Load the simulation parameters.
    """
    p = self._srf["/"+self._prefixb+str(blockid)+"/simulation_parameters"].attrs
    PP = ParameterProvider.ParameterProvider()

    for key, value in p.iteritems():
        PP[key] = pickle.loads(value)

    # Compute some values on top of the given input parameters
    PP.compute_parameters()

    return PP


def update_parameters(self, parameters, blockid="global"):
    params = self.load_parameters(blockid=blockid)
    self.delete_parameters(blockid=blockid)
    params.update_parameters(parameters)
    self.add_parameters(blockid=blockid)
    self.save_parameters(params, blockid=blockid)
