"""The WaveBlocks Project

Compute the autocorrelations of wavefunctions.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2016 R. Bourquin
@license: Modified BSD License
"""

from numpy import array, product, conjugate, floating
from numpy.fft import fftn, ifftn

from WaveBlocksND import BlockFactory
from WaveBlocksND import WaveFunction
from WaveBlocksND import BasisTransformationWF


def compute_autocorrelation(iom, obsconfig=None, blockid=0, eigentrafo=True):
    """Compute the autocorrelation of a wavefunction timeseries.

    :param iom: An :py:class:`IOManager` instance providing the simulation data.
    :param obsconfig: Configuration parameters describing f.e. the inner product to use.
    :type obsconfig: A :py:class:`ParameterProvider` instance.
                     Value has no effect in this class.
    :param blockid: The data block from which the values are read.
    :type blockid: Integer, Default is ``0``
    :param eigentrafo: Whether to make a transformation into the eigenbasis.
    :type eigentrafo: Boolean, default is ``True``.
    """
    parameters = iom.load_parameters()

    # Number of time steps we saved
    timesteps = iom.load_wavefunction_timegrid(blockid=blockid)
    nrtimesteps = timesteps.shape[0]

    # Construct the grid from the parameters
    grid = BlockFactory().create_grid(parameters)

    # Basis transformator
    if eigentrafo is True:
        # The potential used
        Potential = BlockFactory().create_potential(parameters)
        BT = BasisTransformationWF(Potential)
        BT.set_grid(grid)

    # And two empty wavefunctions
    WFo = WaveFunction(parameters)
    WFo.set_grid(grid)

    WFt = WaveFunction(parameters)
    WFt.set_grid(grid)

    # We want to save norms, thus add a data slot to the data file
    iom.add_autocorrelation(parameters, timeslots=nrtimesteps, blockid=blockid)

    # Preconfigure the
    values = iom.load_wavefunction(timestep=0, blockid=blockid)
    values = [values[j, ...] for j in range(parameters["ncomponents"])]
    WFo.set_values(values)

    # Project wavefunction values to eigenbasis
    if eigentrafo is True:
        BT.transform_to_eigen(WFo)

    # Fourier transform the values
    WFo.set_values([fftn(value) for value in WFo.get_values()])

    # Iterate over all timesteps
    for i, step in enumerate(timesteps):
        print(" Computing autocorrelations of timestep %d" % step)

        # Retrieve simulation data
        values = iom.load_wavefunction(timestep=step, blockid=blockid)
        values = [values[j, ...] for j in range(parameters["ncomponents"])]
        WFt.set_values(values)

        # Project wavefunction values to eigenbasis
        if eigentrafo is True:
            BT.transform_to_eigen(WFt)

        # Fourier transform the values
        WFt.set_values([fftn(value) for value in WFt.get_values()])

        # Compute the prefactor
        T = grid.get_extensions()
        N = grid.get_number_nodes()
        prefactor = product(array(T) / array(N).astype(floating)**2)

        # Compute the autocorrelation
        # TODO: Consider splitting into cases `fft` versus `fftn`
        valueso = WFo.get_values()
        valuest = WFt.get_values()
        acs = [prefactor * ifftn(sum(conjugate(valueso[n]) * valuest[n])) for n in range(parameters["ncomponents"])]

        iom.save_autocorrelation(acs, timestep=step, blockid=blockid)
