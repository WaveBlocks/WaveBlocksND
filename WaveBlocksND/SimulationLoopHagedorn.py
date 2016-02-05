"""The WaveBlocks Project

This file contains the main simulation loop
for the homogeneous Hagedorn propagator.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2015, 2016 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND.SimulationLoop import SimulationLoop
from WaveBlocksND.IOManager import IOManager
from WaveBlocksND.TimeManager import TimeManager
from WaveBlocksND.BlockFactory import BlockFactory
from WaveBlocksND.BasisTransformationHAWP import BasisTransformationHAWP

__all__ = ["SimulationLoopHagedorn"]


class SimulationLoopHagedorn(SimulationLoop):
    r"""This class acts as the main simulation loop. It owns a propagator that
    propagates a set of initial values during a time evolution.
    """

    def __init__(self, parameters, resultsfile):
        r"""Create a new simulation loop instance for a simulation
        using the semiclassical Hagedorn wavepacket based propagation
        method.

        :param parameters: The simulation parameters.
        :type parameters: A :py:class:`ParameterProvider` instance.
        :param resultsfile: Path and filename of the hdf5 output file.
        """
        # Keep a reference to the simulation parameters
        self.parameters = parameters

        # The time propagator instance driving the simulation.
        self.propagator = None

        # A `IOManager` instance for saving simulation results.
        self.IOManager = None

        # The time manager
        self._tm = TimeManager(self.parameters)

        # Set up serialization of simulation data
        self.IOManager = IOManager()
        self.IOManager.create_file(resultsfile)

        # Save the simulation parameters
        self.IOManager.add_parameters()
        self.IOManager.save_parameters(parameters)


    def prepare_simulation(self):
        r"""Set up a Hagedorn propagator for the simulation loop. Set the
        potential and initial values according to the configuration.

        :raise: :py:class:`ValueError` For invalid or missing input data.
        """
        # The potential instance
        potential = BlockFactory().create_potential(self.parameters)

        # Project the initial values to the canonical basis
        BT = BasisTransformationHAWP(potential)

        # Finally create and initialize the propagator instance
        # TODO: Attach the "leading_component to the hawp as codata
        # TODO: Clean up this ugly if tree
        if self.parameters["propagator"] == "magnus_split":
            from .MagnusPropagator import MagnusPropagator
            self.propagator = MagnusPropagator(self.parameters, potential)
        elif self.parameters["propagator"] == "semiclassical":
            from .SemiclassicalPropagator import SemiclassicalPropagator
            self.propagator = SemiclassicalPropagator(self.parameters, potential)
        elif self.parameters["propagator"] == "McL42sc":
            from .McL42scPropagator import McL42scPropagator
            self.propagator = McL42scPropagator(self.parameters, potential)
        elif self.parameters["propagator"] == "Pre764sc":
            from .Pre764scPropagator import Pre764scPropagator
            self.propagator = Pre764scPropagator(self.parameters, potential)
        elif self.parameters["propagator"] == "McL84sc":
            from .McL84scPropagator import McL84scPropagator
            self.propagator = McL84scPropagator(self.parameters, potential)
        elif self.parameters["propagator"] == "hagedorn":
            from .HagedornPropagator import HagedornPropagator
            self.propagator = HagedornPropagator(self.parameters, potential)
        else:
            raise NotImplementedError("Unknown propagator type: " + self.parameters["propagator"])

        # Create suitable wavepackets
        chi = self.parameters["leading_component"]

        for packet_descr in self.parameters["initvals"]:
            packet = BlockFactory().create_wavepacket(packet_descr)
            # Transform to canonical basis
            BT.set_matrix_builder(packet.get_innerproduct())
            BT.transform_to_canonical(packet)
            # And hand over
            self.propagator.add_wavepacket((packet, chi))

        # Add storage for each packet
        npackets = len(self.parameters["initvals"])
        slots = self._tm.compute_number_events()
        key = ("q","p","Q","P","S","adQ")

        for i in range(npackets):
            bid = self.IOManager.create_block()
            self.IOManager.add_wavepacket(self.parameters, timeslots=slots, blockid=bid, key=key)

        # Write some initial values to disk
        for packet in self.propagator.get_wavepackets():
            self.IOManager.save_wavepacket_description(packet.get_description())

        if self._tm.is_event(0):
            for packet in self.propagator.get_wavepackets():
                # Pi
                self.IOManager.save_wavepacket_parameters(packet.get_parameters(key=key), timestep=0, key=key)
                # Basis shapes
                for shape in packet.get_basis_shapes():
                    self.IOManager.save_wavepacket_basisshapes(shape)
                # Coefficients
                self.IOManager.save_wavepacket_coefficients(packet.get_coefficients(), packet.get_basis_shapes(), timestep=0)


    def run_simulation(self):
        r"""Run the simulation loop for a number of time steps.
        """
        # The number of time steps we will perform.
        nsteps = self._tm.compute_number_timesteps()

        # Which parameter data to save.
        key = ("q","p","Q","P","S","adQ")

        # Run the prepropagate step
        self.propagator.pre_propagate()
        # Note: We do not save any data here

        # Run the simulation for a given number of timesteps
        for i in range(1, nsteps+1):
            print(" doing timestep "+str(i))

            self.propagator.propagate()

            # Save some simulation data
            if self._tm.is_event(i):
                # Run the postpropagate step
                self.propagator.post_propagate()

                # TODO: Generalize for arbitrary number of wavepackets
                packets = self.propagator.get_wavepackets()
                assert len(packets) == 1

                for packet in packets:
                    # Pi
                    self.IOManager.save_wavepacket_parameters(packet.get_parameters(key=key), timestep=i, key=key)
                    # Basis shapes (in case they changed!)
                    for shape in packet.get_basis_shapes():
                        self.IOManager.save_wavepacket_basisshapes(shape)
                    # Coefficients
                    self.IOManager.save_wavepacket_coefficients(packet.get_coefficients(), packet.get_basis_shapes(), timestep=i)

                # Run the prepropagate step
                self.propagator.pre_propagate()

        # Run the postpropagate step
        self.propagator.post_propagate()
        # Note: We do not save any data here


    def end_simulation(self):
        r"""Do the necessary cleanup after a simulation. For example request the
        :py:class:`IOManager` to write the data and close the output files.
        """
        self.IOManager.finalize()

        print("**************************************************")
        print("***            Simulation Finished             ***")
        print("**************************************************")
