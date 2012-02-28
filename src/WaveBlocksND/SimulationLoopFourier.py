"""The WaveBlocks Project

This file contains the main simulation loop
for the Fourier propagator.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import numpy as np
import scipy as sp


from TensorProductGrid import TensorProductGrid
from PotentialFactory import PotentialFactory
from WaveFunction import WaveFunction
from Initializer import Initializer
from FourierPropagator import FourierPropagator
from SimulationLoop import SimulationLoop
#from IOManager import IOManager

__all__ = ["SimulationLoopFourier"]


class SimulationLoopFourier(SimulationLoop):
    """This class acts as the main simulation loop. It owns a propagator that
    propagates a set of initial values during a time evolution.
    """

    def __init__(self, parameters):
        """Create a new simulation loop instance for a simulation
        using the Fourier propagation method.

        :param parameters: The simulation parameters.
        :type parameters: A :py:class:`ParameterProvider` instance.
        """
        # Keep a reference to the simulation parameters
        self.parameters = parameters

        # The time propagator instance driving the simulation.
        self.propagator = None

        # An `IOManager` instance for saving simulation results.
        self.IOManager = None

        # Which data do we want to save
        self._tm = self.parameters.get_timemanager()

        # TODO: Set up serializing of simulation data
        #self.IOManager = IOManager()
        #self.IOManager.create_file(self.parameters)
        #self.IOManager.create_block()


    def prepare_simulation(self):
        """Set up a Fourier propagator for the simulation loop. Set the
        potential and initial values according to the configuration.
        :raise ValueError: For invalid or missing input data.
        """
        # The potential instance
        potential = PotentialFactory().create_potential(self.parameters)

        # Compute the position space grid points
        grid = TensorProductGrid(self.parameters)

        # Construct initial values
        I = Initializer(self.parameters)
        initialvalues = I.initialize_for_fourier(grid)

        # Project the initial values to the canonical basis
        # TODO: Make basis transformation procedure working with WF objects
        #initialvalues = potential.project_to_canonical(nodes, initialvalues.get_values())

        # Store the initial values in a WaveFunction object
        #IV = WaveFunction(self.parameters)
        #IV.set_grid(grid)
        #IV.set_values(initialvalues)
        IV = initialvalues

        # Finally create and initialize the propagator instace
        self.propagator = FourierPropagator(potential, IV, self.parameters)

        #slots = self._tm.compute_number_saves()

        #self.IOManager.add_grid(self.parameters, blockid="global")
        #self.IOManager.add_fourieroperators(self.parameters)
        #self.IOManager.add_wavefunction(self.parameters, timeslots=slots)

        # Write some initial values to disk
        #self.IOManager.save_grid(nodes, blockid="global")
        #self.IOManager.save_fourieroperators(self.propagator.get_operators())
        #self.IOManager.save_wavefunction(IV.get_values(), timestep=0)


    def run_simulation(self):
        """Run the simulation loop for a number of time steps. The number of steps
        is calculated in the :py:method:`prepare` function."""
        # The number of time steps we will perform.
        nsteps = self._tm.compute_number_timesteps()

        # Run the simulation for a given number of timesteps
        for i in xrange(1, nsteps+1):
            print(" doing timestep "+str(i))

            self.propagator.propagate()

            # Save some simulation data
            # TODO
            #if tm.must_save(i):
            #    self.IOManager.save_wavefunction(self.propagator.get_wavefunction().get_values(), timestep=i)


    def end_simulation(self):
        """Do the necessary cleanup after a simulation. For example request the
        :py:class:`IOManager` to write the data and close the output files.
        """
        # TODO
        #self.IOManager.finalize()
        pass
