import argparse
import os
from numpy import exp, sqrt, linspace

from WaveBlocksND import ParameterLoader
from WaveBlocksND import TimeManager
from WaveBlocksND import HyperCubicShape
from WaveBlocksND import MorseEigenstate
from WaveBlocksND import IOManager
from WaveBlocksND import GlobalDefaults


parser = argparse.ArgumentParser()

parser.add_argument("parametersfile",
                    type = str,
                    help = "The simulation configuration parameters file.")

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
    raise IOError("The results path does not exist: {}".format(args.resultspath))

outputfile = os.path.abspath(os.path.join(args.resultspath, args.outputfile))
parametersfile = os.path.abspath(args.parametersfile)

# Set up the parameter provider singleton
PP = ParameterLoader().load_from_file(args.parametersfile)

# Initial value
MWP = MorseEigenstate(PP['eps'], PP['beta'], PP['V0'])
Kmax = min(PP['Kmax'], MWP.get_number_eigenstates())
K = HyperCubicShape([Kmax])
MWP.set_basis_shapes(K)

C = MWP.get_coefficients()
C[0] = 1.0 / sqrt(2)
C[1] = 1.0 / sqrt(2)
MWP.set_coefficients(C)

TM = TimeManager(PP)

# Set up serialization of simulation data
IOM = IOManager()
IOM.create_file(outputfile)
IOM.create_block()

# Save the simulation parameters
IOM.add_parameters()
IOM.save_parameters(PP)

# Save grid and wavefunction values
IOM.add_grid(PP, blockid="global")
IOM.add_wavefunction(PP, timeslots=None)

# Evaluate initial value
X = linspace(PP['limits'][0][0], PP['limits'][0][1], PP['number_nodes'][0]).reshape(1, -1)
WF = MWP.evaluate_at(X)

IOM.save_grid(X, blockid="global")
if TM.is_event(0):
    IOM.save_wavefunction(WF, timestep=0)

# Propagator
E = MWP.energy_levels()
P = exp(-1.0j * PP['dt'] * E / PP['eps']**2).reshape(-1, 1)

# The number of time steps we will perform
nsteps = TM.compute_number_timesteps()

# Run the simulation for a given number of timesteps
for i in range(1, nsteps + 1):
    print(" doing timestep {}".format(i))

    Cold = MWP.get_coefficients()
    Cnew = P * Cold
    Cnew = MWP.set_coefficients(Cnew)

    if TM.is_event(i):
        WF = MWP.evaluate_at(X)
        IOM.save_wavefunction(WF, timestep=i)
