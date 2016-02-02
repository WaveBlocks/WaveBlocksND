#!/usr/bin/env python
"""The WaveBlocks Project

This file is main script for running simulations with WaveBlocks.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os

from WaveBlocksND import ParameterLoader
from WaveBlocksND import GlobalDefaults

parser = argparse.ArgumentParser()

parser.add_argument("parametersfile",
                    type = str,
                    help = "The simulation configuration parameters file")

parser.add_argument("-o", "--outputfile",
                    type = str,
                    help = "The data file to write the transformed data.",
                    default = GlobalDefaults.file_resultdatafile)

parser.add_argument("-r", "--resultspath",
                    type = str,
                    help = "Path where to put the results.",
                    nargs = "?",
                    default = '.')


args = parser.parse_args()


# Check if the results path exists and assemble output file path
resultspath = os.path.abspath(args.resultspath)

if not os.path.exists(resultspath):
    raise IOError("The results path does not exist: " + args.resultspath)

outputfile = os.path.abspath(os.path.join(args.resultspath, args.outputfile))
parametersfile = os.path.abspath(args.parametersfile)

print("Using configuration from file: " + parametersfile)
print("Storing simulation results at: " + resultspath)
print("Output data file is          : " + outputfile)


# Set up the parameter provider singleton
PA = ParameterLoader().load_from_file(args.parametersfile)

# Print the parameters that apply for this simulation
print(PA)

# Decide which simulation loop to use
if PA["algorithm"] == "fourier":
    from WaveBlocksND import SimulationLoopFourier
    SL = SimulationLoopFourier(PA, resultsfile=outputfile)

elif PA["algorithm"] == "hagedorn":
    from WaveBlocksND import SimulationLoopHagedorn
    SL = SimulationLoopHagedorn(PA, resultsfile=outputfile)

elif PA["algorithm"] == "hagedorn_inhomog":
    from WaveBlocksND import SimulationLoopHagedornInhomogeneous
    SL = SimulationLoopHagedornInhomogeneous(PA, resultsfile=outputfile)

# NOTE: Add new algorithms here

else:
    raise ValueError("Invalid propagator algorithm.")

# Initialize and run the simulation
SL.prepare_simulation()
SL.run_simulation()

# End the simulation, close output files etc.
SL.end_simulation()
