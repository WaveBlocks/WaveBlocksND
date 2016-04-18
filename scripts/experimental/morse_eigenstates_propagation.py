from numpy import arange, zeros, exp, sqrt, floor, linspace, complexfloating, sum, zeros_like, pi, log
from scipy.special import gamma, eval_genlaguerre

from WaveBlocksND import TimeManager, IOManager, ParameterLoader


class MorseWavepacket:

    def __init__(self, eps, beta, V0):
        self._eps = eps
        
        self._beta = beta
        self._V0 = V0

        self._nu = sqrt(8 * V0 / (beta**2 * eps**4))
      
        
    def get_nu(self):
        return self._nu

    
    def get_max_levels(self):
        Kmax = floor((self._nu - 1) / 2)
        return Kmax

    
    def Nn(self, n):
        return sqrt(self._beta * (self._nu - 2 * n - 1) * gamma(n + 1) / gamma(self._nu - n))


    def mu(self, n, x):
        sn = 1 / 2 * (self._nu - 2 * n - 1)
        z = self._nu * exp(-self._beta * x.reshape(-1))
        mun = self.Nn(n) * exp(-z / 2) * z**sn * eval_genlaguerre(n, 2 * sn, z)
        return mun

    
    def _mu0(self, x):
        # Improved version for small epsilon
        delta = ((self._nu - 1) + 1 / (12 * (self._nu - 1) - 1 / (10 * (self._nu - 1))))
        r2 = sqrt(self._beta * delta / (exp(1) * self._nu))
        r4 = ((self._nu - 1) / (2 * pi)) ** 0.25
        ex = log(delta) - 1 - log(self._nu) + (self._nu - 1) / self._nu * self._beta * x + exp(-self._beta * x)
        return r2 * r4 * exp(-self._nu / 2 * ex)

    
    def evaluate_direct(self, K, x):
        Kmax = K.shape[0]
        g = x.shape[1]
        B = zeros((Kmax, g), dtype=complexfloating)

        for n in K:
            B[n, :] = self.mu(n, x)
        return B


    def evaluate_recursive(self, K, x):
        Kmax = K.shape[0]
        g = x.shape[1]
        B = zeros((Kmax, g), dtype=complexfloating)
        
        # mu0
        B[0, :] = self._mu0(x)

        for n in range(0, Kmax - 1):
            rn = sqrt((n + 1) * (self._nu - n - 1))
            ln = sqrt(n * (self._nu - n))
            sn = 1 / 2 * (self._nu - 2 * n - 1)
            pf1 =  1 / rn * sqrt((sn - 1) / sn      ) * (2 * sn)     / (2 * sn + 1)
            pf2 = ln / rn * sqrt((sn - 1) / (sn + 1)) * (2 * sn - 1) / (2 * sn + 1)
            pf1 *= ((4 * sn**2 - 1) * exp(self._beta * x) / self._nu - self._nu)
            # Note: Index wraps around for n = 0
            B[n + 1, :] = pf1 * B[n, :] - pf2 * B[n - 1, :]
            
        return B


    def En(self, n):
        return -1 / 2 * (sqrt(2 * self._V0) - self._eps**2 * abs(self._beta) * (n + 1 / 2))**2





if __name__ == '__main__':
    import argparse
    import os
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
    MWP = MorseWavepacket(PP['eps'], PP['beta'], PP['V0'])
    K = arange(PP['Kmax'])
    C = zeros_like(K, dtype=complexfloating).reshape(-1, 1)
    C[0] = 1.0 / sqrt(2)
    C[1] = 1.0 / sqrt(2)
    #C[11] = 1.0 / sqrt(3)

    TM = TimeManager(PP)

    # Set up serialization of simulation data
    IOM = IOManager()
    IOM.create_file(args.outputfile)
    IOM.create_block()

    # Save the simulation parameters
    IOM.add_parameters()
    IOM.save_parameters(PP)

    # Save grid and wavefunction values
    IOM.add_grid(PP, blockid="global")
    IOM.add_wavefunction(PP, timeslots=None)

    # Evaluate initial value
    X = linspace(PP['limits'][0][0], PP['limits'][0][1], PP['number_nodes'][0]).reshape(1, -1)
    PSI = MWP.evaluate_recursive(K, X)
    WF = sum(C * PSI, axis=0).reshape(1, -1)

    IOM.save_grid(X, blockid="global")
    if TM.is_event(0):
        IOM.save_wavefunction(WF, timestep=0)

    # Propagator
    E = MWP.En(K)
    P = exp(-1.0j * PP['dt'] * E / PP['eps']**2).reshape(-1, 1)

    # The number of time steps we will perform
    nsteps = TM.compute_number_timesteps()

    # Run the simulation for a given number of timesteps
    for i in range(1, nsteps + 1):
        print(" doing timestep {}".format(i))

        PSI = P * PSI
        WF = sum(C * PSI, axis=0).reshape(1, -1)
        IOM.save_wavefunction(WF, timestep=i)
