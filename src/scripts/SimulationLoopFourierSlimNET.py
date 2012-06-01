"""The WaveBlocks Project

This file contains the main simulation loop
for the Fourier propagator.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND import BlockFactory
from WaveBlocksND import Initializer
from WaveBlocksND import BasisTransformationWF
from WaveBlocksND import KineticOperator
from WaveBlocksND import FourierPropagator
from WaveBlocksND import SimulationLoop
from WaveBlocksND import IOManager

__all__ = ["SimulationLoopFourierSlimNET"]


class SimulationLoopFourierSlimNET(SimulationLoop):
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

        # Set up serialization of simulation data
        self.IOManager = IOManager()
        self.IOManager.create_file(self.parameters)
        self.IOManager.create_block()


    def prepare_simulation(self):
        r"""Set up a Fourier propagator for the simulation loop. Set the
        potential and initial values according to the configuration.

        :raise ValueError: For invalid or missing input data.
        """
        # The potential instance
        potential = BlockFactory().create_potential(self.parameters)

        # Compute the position space grid points
        grid = BlockFactory().create_grid(self.parameters)

        # Construct initial values
        I = Initializer(self.parameters)
        initialvalues = I.initialize_for_fourier(grid)

        # Transform the initial values to the canonical basis
        self._BT = BasisTransformationWF(potential)
        self._BT.set_grid(grid)
        self._BT.transform_to_canonical(initialvalues)

        # Finally create and initialize the propagator instace
        self.propagator = FourierPropagator(potential, initialvalues, self.parameters)

        # The operators for computing energies
        KO = KineticOperator(grid)
        KO.calculate_operator(self.parameters["eps"])
        self._opT = KO
        # There is NO eigen transformation involved here
        self._opV = potential.evaluate_eigenvalues_at(grid, as_matrix=True)

        # Write some initial values to disk
        slots = self._tm.compute_number_saves()

        # What data to save?
        # We want to save norms, thus add a data slot to the data file
        self.IOManager.add_norm(self.parameters, timeslots=slots)
        # We want to save norms, thus add a data slot to the data file
        self.IOManager.add_energy(self.parameters, timeslots=slots)

        # Save the norm at the beginning of our simulation
        WF = self.propagator.get_wavefunction()
        norms = WF.norm()
        self.IOManager.save_norm(norms, timestep=0)

        # Save the energies at the beginning of our simulation
        ekinlist = WF.kinetic_energy(self._opT, summed=False)
        epotlist = WF.potential_energy(self._opV, summed=False)
        self.IOManager.save_energy((ekinlist, epotlist), timestep=0)


    def run_simulation(self):
        r"""Run the simulation loop for a number of time steps.
        """
        # The number of time steps we will perform.
        nsteps = self._tm.compute_number_timesteps()

        # Run the simulation for a given number of timesteps
        for i in xrange(1, nsteps+1):
            print(" doing timestep "+str(i))

            self.propagator.propagate()

            # Save some simulation data
            if self._tm.must_save(i):
                # Calculate the norm of the wave functions
                WF = self.propagator.get_wavefunction()
                norms = WF.norm()
                self.IOManager.save_norm(norms, timestep=i)
                # Compute and save enrgies
                ekinlist = WF.kinetic_energy(self._opT, summed=False)
                epotlist = WF.potential_energy(self._opV, summed=False)
                self.IOManager.save_energy((ekinlist, epotlist), timestep=i)


    def end_simulation(self):
        """Do the necessary cleanup after a simulation. For example request the
        :py:class:`IOManager` to write the data and close the output files.
        """
        self.IOManager.finalize()
