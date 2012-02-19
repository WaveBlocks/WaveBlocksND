"""The WaveBlocks Project

This file contains the Fourier propagator class. The wavefunction
:math:`Psi` is propagated in time with a strang splitting of the
exponential :math:`\exp(-\frac{i}{\varepsilon^2} \tau H)`.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import arange, append, exp, zeros, complexfloating, squeeze
from numpy.fft import fftn, ifftn

from Propagator import Propagator
from KineticOperator import KineticOperator


class FourierPropagator(Propagator):
    # """This class can numerically propagate given initial values $\Ket{\Psi}$ in
    # a potential surface $V\ofs{x}$. The propagation is done with a Strang splitting
    # of the time propagation operator."""

    def __init__(self, grid, potential, initial_values, para):
        # """Initialize a new I{FourierPropagator} instance. Precalculate also the
        # grid and the propagation operators.
        # @param potential: The potential the state $\Ket{\Psi}$ feels during the time propagation.
        # @param initial_values: The initial values $\Ket{\Psi\ofs{t=0}}$ given in the canonical basis.
        # @raise ValueError: If the number of components of $\Ket{\Psi}$ does not
        # match the number of energy levels $\lambda_i$ of the potential.
        # """
        #: The embedded I{MatrixPotential} instance representing the potential $V$.
        self.potential = potential

        #: The initial values of the components $\psi_i$ sampled at the given nodes.
        self.Psi = initial_values

        if self.potential.get_number_components() != self.Psi.get_number_components():
            raise ValueError("Potential dimension and number of states do not match")

        #: The position space nodes $\gamma$.
        self.grid = grid

        #: The potential operator $V$ defined in position space.
        self.V = self.potential.evaluate_at(grid)

        self.KO = KineticOperator(grid)
        self.KO.calculate_operator(para["dt"], para["eps"])

        #: The kinetic operator $T$ defined in momentum space.
        self.T = self.KO.evaluate_at()

        #: Exponential $\exp\ofs{T}$ of $T$ used in the Strang splitting.
        self.TE = self.KO.evaluate_exponential_at()

        self.potential.calculate_exponential(-0.5j * para["dt"] / para["eps"]**2)
        #: Exponential $\exp\ofs{V}$ of $V$ used in the Strang splitting.
        self.VE = [ self.potential.evaluate_exponential_at(grid) ]


    def __str__(self):
        """Prepare a printable string representing the I{FourierPropagator} instance."""
        return "Fourier propagator for " + str(self.potential.get_number_components()) + " components."


    def get_number_components(self):
        """@return: The number of components of $\Ket{\Psi}$."""
        return self.potential.get_number_components()


    def get_potential(self):
        """@return: The I{MatrixPotential} instance used for time propagation."""
        return self.potential


    def get_wavefunction(self):
        """@return: The I{WaveFunction} instance that stores the current wavefunction data."""
        return self.Psi


    # def get_operators(self):
    #     """@return: Return the numerical expressions of the propagation
    #     operators $T$ and $V$.
    #     """
    #     return (self.T, self.V)


    def propagate(self):
        """Given the wavefunction values :math:`\Psi(\Gamma)` at time :math:`t`, calculate
        new values :math:`\Psi^\prime(\Gamma)` at time :math:`t + \tau`. We perform exactly
        one single timestep of size :math:`\tau` within this function.
        """
        # How many components does Psi have
        N = self.Psi.get_number_components()

        # Unpack the values from the current WaveFunction
        values = self.Psi.get_values()

        tmp = [ zeros(value.shape, dtype=complexfloating) for value in values ]

        # The first step with the potential
        for row in xrange(0, N):
            for col in xrange(0, N):
                tmp[row] = tmp[row] + self.VE[row*N+col] * values[col]

        # Go to Fourier space
        tmp = [ fftn(component) for component in tmp ]

        # Apply the kinetic operator
        tmp = [ self.TE * component for component in tmp ]

        # Go back to real space
        tmp = [ ifftn(component) for component in tmp ]

        # The second step with the potential
        values = [ zeros(component.shape, dtype=complexfloating) for component in tmp ]
        for row in xrange(0, N):
            for col in xrange(0, N):
                values[row] = values[row] + self.VE[row*N+col] * tmp[col]

        # Pack values back to WaveFunction object
        # TODO: Consider squueze(.) of data befor repacking
        self.Psi.set_values(values)
