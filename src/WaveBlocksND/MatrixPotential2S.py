r"""The WaveBlocks Project

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
from Grid import Grid
from GridWrapper import GridWrapper

__all__ = ["MatrixPotential2S"]


class MatrixPotential2S(MatrixPotential):
    r"""This class represents a matrix potential :math:`V(x)`. The potential is
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

        # The the potential, symbolic expressions and evaluatable functions
        assert expression.shape == (2,2)

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

        # The cached Jacobian of the eigenvalues, symbolic expressions and evaluatable functions
        self._jacobian_s = None
        self._jacobian_n = None

        # The cached Hessian of the eigenvalues, symbolic expressions and evaluatable functions
        self._hessian_s = None
        self._hessian_n = None


    def _grid_wrap(self, grid):
        # TODO: Consider additional input types for "nodes":
        #       list of numpy ndarrays, list of single python scalars
        if not isinstance(grid, Grid):
            if not type(grid) is numpy.ndarray:
                grid = numpy.atleast_2d(grid)
            grid = GridWrapper(grid)
        return grid


    def evaluate_at(self, grid, entry=None, as_matrix=True):
        r"""Evaluate the potential :math:`V(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the potential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The indices :math:`(i,j)` of the component :math:`V_{i,j}(x)`
                      we want to evaluate or `None` to evaluate all entries.
        :type entry: A python tuple of two integers.
        :param as_matrix: Dummy parameter which has no effect here.
        :return: A list containing 4 numpy ndarrays of shape :math:`(1, |\Gamma|)`.
        """
        grid = self._grid_wrap(grid)

        # Determine which entries to evaluate
        if entry is not None:
            (row, col) = entry
            entries = [ row * self._number_components + col ]
        else:
            N = self._number_components
            entries = [ row * N + col for row in xrange(N) for col in xrange(N) ]
            #entries = [ n for n in xrange(N*N) ]

        # Evaluate all entries specified
        result = []

        nodes = grid.get_nodes(split=True)

        for index in entries:
            # Evaluate the potential at the given nodes
            values = self._potential_n[index](*nodes)

            # Test for potential beeing constant
            if numpy.atleast_1d(values).shape == (1,):
                values = values * numpy.ones(grid.get_number_nodes(), dtype=numpy.complexfloating)

            # Put the result in correct shape (1, #gridnodes)
            N = grid.get_number_nodes(overall=True)
            result.append( values.reshape((1,N)) )

        # TODO: Consider unpacking single ndarray iff entry != None
        if entry is not None:
            result = result[0]

        return result


    def calculate_eigenvalues(self):
        r"""Calculate the two eigenvalues :math:`\lambda_i(x)` of the potential :math:`V(x)`.
        We can do this by symbolic calculations. The multiplicities are taken into account.
        Note: This function is idempotent and the eigenvalues are memoized for later reuse.
        """
        if self._eigenvalues_s is not None:
            return

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
        r"""Evaluate the eigenvalues :math:`\lambda_i(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid :math:`\Gamma` containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvalues at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The indices :math:`(i,j)` of the component :math:`\Lambda_{i,j}(x)`
                      we want to evaluate or `None` to evaluate all entries. If :math:`j = i`
                      then we evaluate the eigenvalue :math:`\lambda_i(x)`.
        :type entry: A python tuple of two integers.
        :param as_matrix: Whether to include the off-diagonal zero entries of
                          :math:`\Lambda_{i,j}(x)` in the return value.
        :return: A list containing the numpy ``ndarray``, all of shape :math:`(1, |\Gamma|)`.
        """
        if self._eigenvalues_n is None:
            self.calculate_eigenvalues()

        grid = self._grid_wrap(grid)

        # Determine which entries to evaluate
        if entry is not None:
            # Single entry only
            entries = [ entry ]
        else:
            N = self._number_components
            if as_matrix is True:
                # All entries
                entries = [ (row,col) for row in xrange(N) for col in xrange(N) ]
            else:
                # Diagonal entries only
                entries = [ (row,row) for row in xrange(N) ]

        # Compute all diagonal entries first
        diags = [ entry[0] for entry in entries if entry[0] == entry[1] ]

        nodes = grid.get_nodes(split=True)

        tmp = []

        for index in diags:
            # Evaluate the eigenvalue at the given nodes
            values = self._eigenvalues_n[index]( *nodes )

            # Test for eigenvalue beeing constant
            if numpy.atleast_1d(values).shape == (1,):
                values = values * numpy.ones(grid.get_number_nodes(), dtype=numpy.complexfloating)

            tmp.append(values)

        # Sort the eigenvalues pointwise. We can do this because we
        # assume that the different eigenvalues never cross.
        # TODO: Sort will fail iff energy levels really cross!
        N = len(diags)

        if N > 1:
            tmp = numpy.vsplit(numpy.sort(numpy.vstack(tmp), axis=0), N)
            # Take in descending order and reshape
            tmp = [ item.reshape((1,grid.get_number_nodes(overall=True))) for item in reversed(tmp) ]

        # Compose the result for all entries specified
        result = []

        for entry in entries:
            (row, col) = entry
            if row == col:
                result.append( tmp[row] )
            else:
                # Evaluate an off-diagonal entry which equals zero by definition
                result.append( numpy.zeros((1,grid.get_number_nodes(overall=True)), dtype=numpy.complexfloating) )

        # TODO: Consider unpacking single ndarray iff entry != None
        if entry is not None:
            result = result[0]

        return result


    def calculate_eigenvectors(self):
        r"""Calculate the two eigenvectors :math:`\nu_i(x)` of the potential :math:`V(x)`.
        We can do this by symbolic calculations.
        Note: This function is idempotent and the eigenvectors are memoized for later reuse.
        """
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

        # Attention, the components get listed in columns-wise order!
        for vector in self._eigenvectors_s:
            self._eigenvectors_n.append([ sympy.lambdify(self._all_variables, component, "numpy")
                                          for component in vector ])


    def evaluate_eigenvectors_at(self, grid, entry=None):
        r"""Evaluate the two eigenvectors :math:`\nu_i(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvectors at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The index :math:`i` of the eigenvector :math:`\nu_i(x)`
                      we want to evaluate or `None` to evaluate all eigenvectors.
        :type entry: A singly python  integer.
        :return: A list containing the numpy ndarrays, all of shape :math:`(N, |\Gamma|)`.
        """
        # TODO: Rethink about the 'entry' parameter here. Do we need it?
        if self._eigenvectors_n is None:
            self.calculate_eigenvectors()

        grid = self._grid_wrap(grid)

        nodes = grid.get_nodes(split=True)
        # Assure real values as atan2 is only defined for real values!
        nodes = map(numpy.real, nodes)

        # Evaluate the eigenvectors on the grid
        result = []

        for vector in self._eigenvectors_n:
            tmp = numpy.zeros((self._number_components, grid.get_number_nodes(overall=True)), dtype=numpy.complexfloating)
            for index in xrange(self._number_components):
                tmp[index,:] = vector[index](*nodes)
            result.append(tmp)

        return tuple(result)


    def calculate_exponential(self, factor=1):
        r"""Calculate the matrix exponential :math:`\exp(\alpha V)`. In the case
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

        try:
            D = sympy.simplify(D)
            t = sympy.simplify(t)
        except:
            pass

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

        self._exponential_s = M
        self._exponential_n = tuple([ sympy.lambdify(self._all_variables, item, "numpy")
                                      for item in self._exponential_s ])


    def evaluate_exponential_at(self, grid):
        r"""Evaluate the exponential of the potential matrix :math:`V(x)` on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the exponential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :return: The numerical approximation of the matrix exponential at the given grid nodes.
                 A list contains the exponentials for all entries :math:`(i,j)`, each having the
                 same shape as the grid.
        """
        grid = self._grid_wrap(grid)

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


    def calculate_jacobian(self):
        r"""Calculate the Jacobian matrix :math:`\nabla \lambda_i` of the potential's
        eigenvalues :math:`\Lambda(x)` with :math:`x \in \mathbb{R}^D`. For potentials
        which depend only one variable, this equals the first derivative and :math:`D=1`.
        Note that this function is idempotent.
        """
        if self._jacobian_s is None:
            self.calculate_eigenvalues()

            self._jacobian_s = []
            # TODO: Add symbolic simplification
            for ew in self._eigenvalues_s:
                tmp = sympy.Matrix([[ew]])
                self._jacobian_s.append( tmp.jacobian(self._all_variables).T )

            self._jacobian_n = []

            # Attention, the components get listed in columns-wise order!
            for jacobian in self._jacobian_s:
                self._jacobian_n.append([ sympy.lambdify(self._all_variables, component, "numpy")
                                          for component in jacobian ])


    def evaluate_jacobian_at(self, grid, component=None):
        r"""Evaluate the list of Jacobian matrices :math:`\nabla \lambda_i(x)` at some grid
        nodes :math:`\Gamma`.

        :param grid: The grid nodes :math:`\Gamma` the Jacobian gets evaluated at.
        :param component: Dummy parameter that has no effect here.
        :return: The value of the potential's Jacobian at the given nodes. The result
                 is a list of ``ndarray`` each of shape :math:`(D,1)` is we evaluate
                 at a single grid node or of shape :math:`(D,|\Gamma|)`
                 if we evaluate at multiple nodes simultaneously.
        """
        # TODO: Rethink about the 'component' parameter here. Do we need it?
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes(split=True)

        D = self._dimension
        N = grid.get_number_nodes(overall=True)

        jacobians = []

        for jacobian in self._jacobian_n:
            J = numpy.zeros((D,N))
            for index, component in enumerate(jacobian):
                J[index, :] = component(*nodes)
            jacobians.append(J)

        return jacobians


    def calculate_hessian(self):
        r"""Calculate the Hessian matrix :math:`\nabla^2 \lambda_i` of the potential's
        eigenvalues :math:`\Lambda(x)` with :math:`x \in \mathbb{R}^D`. For potentials
        which depend only one variable, this equals the second derivative and :math:`D=1`.
        Note that this function is idempotent.
        """
        if self._hessian_s is None:
            self.calculate_eigenvalues()

            self._hessian_s = []
            # TODO: Add symbolic simplification
            for ew in self._eigenvalues_s:
                self._hessian_s.append( sympy.hessian(ew, self._all_variables) )

            self._hessian_n = []

            for hessian in self._hessian_s:
                self._hessian_n.append([ sympy.lambdify(self._all_variables, entry, "numpy")
                                         for entry in hessian ])


    def evaluate_hessian_at(self, grid, component=None):
        r"""Evaluate the list of Hessian matrices :math:`\nabla^2 \lambda_i(x)` at some grid
        nodes :math:`\Gamma`.

        :param grid: The grid nodes :math:`\Gamma` the Hessian gets evaluated at.
        :param component: Dummy parameter that has no effect here.
        :return: The value of the potential's Hessian at the given nodes. The result
                 is an ``ndarray`` of shape :math:`(D,D)` is we evaluate at a single
                 grid node or of shape :math:`(|\Gamma|,D,D)` if we evaluate at multiple
                 nodes simultaneously.
        """
        # TODO: Rethink about the 'component' parameter here. Do we need it?
        grid = self._grid_wrap(grid)
        nodes = grid.get_nodes(split=True)

        D = self._dimension
        N = grid.get_number_nodes(overall=True)

        hessians = []

        for hessian in self._hessian_n:
            H = numpy.zeros((N,D,D))

            for row in xrange(D):
                for col in xrange(D):
                    H[:, row, col] = hessian[row*D+col](*nodes)

            # 'squeeze' would be dangerous here, make sure it works in the 1D case too
            if N == 1:
                H = H[0,:,:]

            hessians.append(H)

        return hessians
