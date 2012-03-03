"""The WaveBlocks Project

This file contains code for the representation of potentials :math:`V(x)`
that contain exactly two energy levels :math:`\lambda_i`. The number
of space dimensions can be arbitrary, :math:`x \in \mathbb{R}^D`.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import sympy
import numpy

from MatrixPotential import MatrixPotential

__all__ = ["MatrixPotential2S"]


class MatrixPotential2S(MatrixPotential):
    """This class represents a matrix potential :math:`V(x)`. The potential is
    given as an analytic :math:`2 \times 2` matrix expression. Some symbolic
    calculations with the potential are supported.
    """

    def __init__(self, expression, variables):
        """Create a new :py:class:`MatrixPotential2S` instance for a given
        potential matrix :math:`V(x)`.

        :param expression:The mathematical expression representing the potential.
        :type expressiom: A `Sympy` matrix type.
        :param variables: The variables corresponding to the space dimensions.
        :type variables: A list of `Sympy` symbols.
        """
        # This class handles potentials with two energy levels.
        self._number_components = 2

        # The variables that represents position space. The order matters!
        self._all_variables = variables

        # The dimension of position space.
        self._dimension = len(variables)

        # Find active variables and constants
        self._variables = []
        self._constants = []
        self._isconstant = []

        for entry in expression:
            # Determine the variables and constants for each component
            local_vars = [ v for v in entry.atoms(sympy.Symbol) ]
            local_consts = [ c for c in variables if not c in local_vars ]
            self._variables.append(local_vars)
            self._constants.append(local_consts)

            # Does the potential depend on all space variables?
            # True:  V_i,j is constant along all axes
            # None:  V_i,j is constant along some axes
            # False: V_i,j is not constant
            if len(local_consts) == 0:
                self._isconstant.append(False)
            elif len(local_vars) == 0:
                self._isconstant.append(True)
            else:
                self._isconstant.append(None)

        # The the potential, symbolic expressions and evaluatable functions
        self._potential_s = expression
        self._potential_n = tuple([ sympy.lambdify(self._all_variables, entry, "numpy") for entry in self._potential_s ])

        # The cached eigenvalues, symbolic expressions and evaluatable functions
        self._eigenvalues_s = None
        self._eigenvalues_n = None

        # The cached eigenvectors, symbolic expressions and evaluatable functions
        self._eigenvectors_s = None
        self._eigenvectors_n = None

        # The cached exponential, symbolic expressions and evaluatable functions
        self._exponential_s = None
        self._exponential_n = None


    def evaluate_at(self, grid, entry=None, as_matrix=True):
        """Evaluate the potential :math:`V(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the potential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The indices :math:`(i,j)` of the component :math:`V_{i,j}(x)`
                      we want to evaluate or `None` to evaluate all entries.
        :type entry: A python tuple of two integers.
        :param as_matrix: Dummy parameter which has no effect here.
        :return: A list containing 4 numpy ndarrays of shape :math:`(N_1, ... ,N_D)`.
        """
        # TODO: Consider additional input types for "grid":
        #       list of numpy ndarrays, list of single python scalars
        #       -> [ atleast_1d(i) for i in values ]
        # if isinstace(grid, Grid):
        #    ...
        # else:
        #    ...

        # Determine which entries to evaluate
        if entry is not None:
            (row, col) = entry
            items = [ row * self._number_components + col ]
        else:
            N = self._number_components
            items = [ row * N + col for row in xrange(N) for col in xrange(N) ]

        # Evaluate all entries specified
        result = []

        for item in items:
            if self._isconstant[item] is False:
                # We can use numpy broadcasting
                result.append( self._potential_n[item](*grid.get_axes()) )
            elif self._isconstant[item] is None:
                values = self._potential_n[item](*grid.get_nodes(split=True))
                # Make sure we work with (N_1, ..., N_D) shaped arrays
                result.append( values.reshape(grid.get_number_nodes()) )
            elif self._isconstant[item] is True:
                values = self._potential_n[item](*grid.get_axes())
                result.append( values * numpy.ones(grid.get_number_nodes(), dtype=numpy.floating) )

        # TODO: Consider unpacking single ndarray iff entry != None
        return result


    def calculate_eigenvalues(self):
        """Calculate the two eigenvalues :math:`\lambda_i(x)` of the potential :math:`V(x)`.
        We can do this by symbolic calculations. The multiplicities are taken into account.
        Note: This function is idempotent and the eigenvalues are memoized for later reuse.
        """
        # TODO: Consider using numerical techniques!

        # Symbolic formula for the eigenvalues of a general 2x2 matrix
        T = self._potential_s.trace()
        D = self._potential_s.det()

        l1 = (T + sympy.sqrt(T**2 - 4*D)) * sympy.Rational(1,2)
        l2 = (T - sympy.sqrt(T**2 - 4*D)) * sympy.Rational(1,2)

        # Symbolic simplification may fail
        try:
            l1 = sympy.simplify(l1)
            l2 = sympy.simplify(l2)
        except:
            pass

        # The symbolic expressions for the eigenvalues
        self._eigenvalues_s = (l1, l2)

        # The numerical functions for the eigenvalues
        self._eigenvalues_n = tuple([ sympy.lambdify(self._all_variables, item, "numpy") for item in self._eigenvalues_s ])


    def evaluate_eigenvalues_at(self, grid, entry=None, as_matrix=False):
        """Evaluate the eigenvalues :math:`\lambda_i(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvalues at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The indices :math:`(i,j)` of the component :math:`\Lambda_{i,j}(x)`
                      we want to evaluate or `None` to evaluate all entries. If :math:`j = i`
                      then we evaluate the eigenvalue :math:`\lambda_i(x)`.
        :type entry: A python tuple of two integers.
        :param as_matrix: Whether to include the off-diagonal zero entries of
                          :math:`\Lambda_{i,j}(x)` in the return value.
        :return: A list containing the numpy ndarrays, all of shape :math:`(N_1, ... ,N_D)`.
        """
        if self._eigenvalues_n is None:
            self.calculate_eigenvalues()

        # Evaluate the eigenvalues on the grid nodes
        # The complex numbers are to avoid issues with sqrt etc returning nans
        nodes = grid.get_nodes(split=True)

        # TODO: Fix for constant eigenvalues!
        tmp = [ numpy.real(f(*nodes)) for f in self._eigenvalues_n ]

        # Sort the eigenvalues pointwise. We can do this because we
        # assume that the different eigenvalues never cross.
        tmp = numpy.sort(numpy.vstack(tmp), axis=0)
        tmp = [ tmp[i,:] for i in reversed(xrange(self._number_components)) ]

        # Make sure we work with (N_1, ..., N_D) shaped arrays
        tmp = [ t.reshape(grid.get_number_nodes()) for t in tmp ]

        # TODO: Make this more efficient, do not evaluate unused components
        if not entry is None:
            (row, col) = entry
            if row == col:
                result = tmp[row]
            else:
                result = numpy.zeros(grid.get_number_nodes(), dtype=numpy.floating)
        elif as_matrix is True:
            result = []
            for row in xrange(self._number_components):
                for col in xrange(self._number_components):
                    if row == col:
                        result.append(tmp[row])
                    else:
                        result.append( numpy.zeros(grid.get_number_nodes(), dtype=numpy.floating) )
        else:
            result = tmp

        return result


    def calculate_eigenvectors(self):
        """Calculate the two eigenvectors :math:`\nu_i(x)` of the potential :math:`V(x)`.
        We can do this by symbolic calculations.
        Note: This function is idempotent and the eigenvectors are memoized for later reuse.
        """
        # TODO: Consider using pure numerical techniques

        # Assumption: The matrix is symmetric
        # TODO: Consider generalization for arbitrary 2x2 matrices?
        V1 = self._potential_s[0,0]
        V2 = self._potential_s[0,1]

        theta = sympy.Rational(1,2) * sympy.atan2(V2,V1)

        # The two eigenvectors
        upper = sympy.Matrix([[ sympy.cos(theta)], [sympy.sin(theta)]])
        lower = sympy.Matrix([[-sympy.sin(theta)], [sympy.cos(theta)]])

        # The symbolic expressions for the eigenvectors
        self._eigenvectors_s = (upper, lower)

        # The numerical functions for the eigenvectors
        self._eigenvectors_n = []

        for vector in self._eigenvectors_s:
            self._eigenvectors_n.append([ sympy.lambdify(self._all_variables, component, "numpy") for component in vector ])


    def evaluate_eigenvectors_at(self, grid, entry=None):
        """Evaluate the two eigenvectors :math:`\nu_i(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvectors at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The index :math:`i` of the eigenvector :math:`\nu_i(x)`
                      we want to evaluate or `None` to evaluate all eigenvectors.
        :type entry: A singly python  integer.
        :return: A list containing the numpy ndarrays, all of shape :math:`(N_1, ... ,N_D, 1)`.
        """
        if self._eigenvectors_n is None:
            self.calculate_eigenvectors()

        # Retrieve the grid nodes
        nodes = grid.get_nodes(split=True)
        nodes = [ numpy.array(n, dtype=numpy.complexfloating) for n in nodes ]
        # Assure real values as atan2 is only defined for real values!
        nodes = map(numpy.real, nodes)

        # Evaluate the eigenvectors on the grid
        result = []

        for vector in self._eigenvectors_n:
            tmp = numpy.zeros((self._number_components, grid.get_number_nodes(overall=True)), dtype=numpy.floating)
            for index in xrange(self._number_components):

                tmp[index,:] = vector[index](*nodes)
            result.append(tmp)

        return tuple(result)


    def calculate_exponential(self, factor=1):
        """Calculate the matrix exponential :math:`\exp(\alpha V)`. In the case
        of this class the matrix is of size :math:`2 \times 2` thus the exponential
        can be calculated analytically for a general matrix.
        Note: This function is idempotent.

        :param factor: The prefactor :math:`\alpha` in the exponential.
        """
        M = factor * self._potential_s
        a = M[0,0]
        b = M[0,1]
        c = M[1,0]
        d = M[1,1]

        D = sympy.sqrt((a-d)**2 + 4*b*c)/2
        t = sympy.exp((a+d)/2)

        M = sympy.Matrix([[0,0],[0,0]])

        if sympy.Eq(D,0):
            # special case
            M[0,0] = t * (1 + (a-d)/2)
            M[0,1] = t * b
            M[1,0] = t * c
            M[1,1] = t * (1 - (a-d)/2)
        else:
            # general case
            M[0,0] = t * (sympy.cosh(D) + (a-d)/2 * sympy.sinh(D)/D)
            M[0,1] = t * (b * sympy.sinh(D)/D)
            M[1,0] = t * (c * sympy.sinh(D)/D)
            M[1,1] = t * (sympy.cosh(D) - (a-d)/2 * sympy.sinh(D)/D)

        print(M)

        self._exponential_s = M
        self._exponential_n = tuple([ sympy.lambdify(self._all_variables, item, "numpy") for item in self._exponential_s ])


    def evaluate_exponential_at(self, grid):
        """Evaluate the exponential of the potential matrix :math:`V(x)` on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the exponential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :return: The numerical approximation of the matrix exponential at the given grid nodes.
                 A list contains the exponentials for all entries :math:`(i,j)`, each having the
                 same shape as the grid.
        """
        tmp = [ f(*grid.get_nodes(split=True)) for f in self._exponential_n ]
        tmp = [ numpy.array(n, dtype=numpy.complexfloating) for n in tmp ]

        # TODO: Better fix for globally constant functions
        tmp2 = []

        for t in tmp:
            if numpy.array(t).ndim == 0:
                tmp2.append( t * numpy.ones(grid.get_number_nodes()) )
            else:
                tmp2.append( t.reshape(grid.get_number_nodes()) )

        return tuple(tmp2)
