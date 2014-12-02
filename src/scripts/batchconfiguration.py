"""The WaveBlocks Project

Default configuration of which scripts are run in the batch loop.
Change the content of the lists as you like but never rename the
variables.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2014 R. Bourquin
@license: Modified BSD License
"""

# The scripts in this list are called once before the simulations are run.
call_prerun = []

# All scripts in this list are called for each simulation configuration
# and with the configuration file as first command line argument
call_simulation = ["Main.py"]

# All scripts in this list are called for each simulation configuration
# but without additional arguments. They can assume that the simulation
# results data file is available at the standard location (default: 'simulation_results.hdf5')
call_for_each = [ # Some useful defaults:
                 "ComputeNormsNET.py",
                 "ComputeEnergiesNET.py",
                 "ComputeAutocorrelationNET.py",

                 "PlotNorms.py",
                 "PlotEnergies.py",
                 "PlotAutocorrelations.py",
                 "PlotWavepacketParameters1D.py",
                 #"PlotWavepacketParametersDD.py",
                 ]

# The scripts in this list are called once after all simulations are finished
# and the results were moved to the final location (default './results/*').
# Put all scripts that do comparisons between different simulations in here.
call_postrun = []

# A complete list of output files that are generated for each simulation run.
# All files matching the patterns in the following list are copied to the
# directory containing the results (default './results/*'). This avoids
# overwriting during the next simulation run.
result_files = ["simulation_results.hdf5",
                "simulation_results_eigen.hdf5",
                "*.pdf",
                "*.png"]
