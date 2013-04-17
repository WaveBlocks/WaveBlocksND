"""The WaveBlocks Project


@author: R. Bourquin
@copyright: Copyright (C) 2013 R. Bourquin
@license: Modified BSD License
"""

from numpy import (zeros, ones, sum, diag, squeeze,  conjugate, transpose, dot,
                   zeros_like, ones_like, product, indices, flipud, array)
from scipy import exp, sqrt
from scipy.linalg import inv, schur

from Quadrature import Quadrature

__all__ = ["NSDInhomogeneous"]


class NSDInhomogeneous(Quadrature):
    r"""
    """

    def __init__(self, QR=None):
        # Pure convenience to allow setting of quadrature rule in constructor
        if QR is not None:
            self.set_qr(QR)
        else:
            self._QR = None


    def __str__(self):
        return "Inhomogeneous numerical steepest descent using a " + str(self._QR)


    def get_description(self):
        r"""Return a description of this quadrature object.
        A description is a ``dict`` containing all key-value pairs
        necessary to reconstruct the current instance. A description
        never contains any data.
        """
        d = {}
        d["type"] = "NSDInhomogeneous"
        d["qr"] = self._QR.get_description()
        return d


    def initialize_packet(self, pacbra, packet=None):
        r"""Provide the wavepacket parts of the inner product to evaluate.
        Since the quadrature is inhomogeneous different wavepackets can be
        used for the 'bra' as well as the 'ket' part.

        :param pacbra: The packet that is used for the 'bra' part.
        :param packet: The packet that is used for the 'ket' part.
        """
        # Allow to ommit the ket if it is the same as the bra
        if packet is None:
            packet = pacbra

        self._pacbra = pacbra
        self._packet = packet

        # TODO: Handle non-oscillatory polynomial part
        self._envelope = lambda V: ones_like(V[0,:])


    def initialize_operator(self, operator=None):
        r"""Provide the operator part of the inner product to evaluate.
        This function initializes the operator used for quadratures.
        For nasty technical reasons there are two functions for
        setting up the operators.

        :param operator: The operator of the inner product.
                         If 'None' a suitable identity is used.
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        #       For this, 'operator' must support the 'component=(r,c)' option.
        # Operator is None is interpreted as identity transformation
        if operator is None:
            self._operator = lambda nodes, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
        else:
            self._operator = operator


    def initialize_operator_matrix(self, operator=None):
        r"""
        Provide the operator part of the inner product to evaluate.
        This function initializes the operator used for building matrices.
        For nasty technical reasons there are two functions for
        setting up the operators.

        :param operator: The operator of the inner product.
                         If 'None' a suitable identity is used.
        """
        # TODO: Make this more efficient, only compute values needed at each (r,c) step.
        # For this, 'operator' must support the 'entry=(r,c)' option.
        N  = self._packet.get_number_components()
        if operator is None:
            # Operator is None is interpreted as identity transformation
            self._operatorm = lambda nodes, entry=None: ones((1,nodes.shape[1])) if entry[0] == entry[1] else zeros((1,nodes.shape[1]))
        else:
            self._operatorm = operator


    def build_bilinear(self, Pibra, Piket):
        r"""Convert the oscillator :math:`-\overline{g_k} + g_l`
        occuring in the integral :math:`\langle\phi_k\left[\Pi_k\right] | \phi_l\left[\Pi_l\right]\rangle`
        into a bilinear form
        :math:`g(x) = \underline{x}^{\mathrm{H}} \mathbf{A} \underline{x} + \underline{b}^{\mathrm{T}} \underline{x} + c`.

        :param Pibra: The parameters :math:`\Pi_k = (\underline{q_k}, \underline{p_k}, \mathbf{Q_k}, \mathbf{P_k})` from the 'bra' packet.
        :param Piket: The parameters :math:`\Pi_l = (\underline{q_l}, \underline{p_l}, \mathbf{Q_l}, \mathbf{P_l})` from the 'ket' packet.
        :return: Three arrays: a matrix :math:`\mathbf{A}` of shape :math:`D \times D`,
                 a vector :math:`\underline{b}` of shape :math:`D \times 1` and a scalar value :math:`c`.
        """
        qk, pk, Qk, Pk = Pibra
        ql, pl, Ql, Pl = Piket

        Gk = dot(Pk, inv(Qk))
        Gl = dot(Pl, inv(Ql))

        # Merge into a single oscillator
        A = 0.5 * (Gl - conjugate(transpose(Gk)))
        b = (0.5 * (  dot(Gk, qk)
                    - dot(conjugate(transpose(Gl)), ql)
                    + dot(transpose(Gk), conjugate(qk))
                    - dot(conjugate(Gl), conjugate(ql))
                   )
             + (pl - conjugate(pk))
            )
        b = conjugate(b)
        c = (0.5 * (  dot(conjugate(transpose(ql)), dot(Gl, ql))
                    - dot(conjugate(transpose(qk)), dot(conjugate(transpose(Gk), qk))))
                 + (dot(conjugate(transpose(qk)),pk) - dot(conjugate(transpose(pl)),ql))
            )
        return A, b, c


    def update_oscillator(T, i):
        r"""Given a upper triangular matrix :math:`\mathbf{T} \in \mathbb{C}^{D \times D}` representing
        the oscillator :math:`g(x) = \underline{x}^{\mathrm{H}} \mathbf{T} \underline{x}`, update
        its entries according to:

        .. math::

            t_{j,j} & := t_{j,j} - \frac{t_{i,j}^2}{4 t_{i,i}} \\
            t_{j,k} & := t_{j,k} - \frac{t_{i,j} t_{i,k}}{2 t_{i,i}} \,, \quad k > j \,.

        :param: The matrix :math:`\mathbf{T}` we want to update. (Note that the matrix
                is not modified and a copy returned.)
        :param: The row :math:`0 \leq i \leq D-1` which is taken as base for the updates.
        :return: The updated matrix :math:`\mathbf{T}`.
        """
        Ti = T.copy()
        rr, cc = T.shape

        if T[i-1,i-1] == 0:
            print("Warning: 'update_oscillator' encountered a RESIDUE situation!")
            return Ti

        # Diagonal Elements
        for j in xrange(i,rr):
            Ti[j,j] = T[j,j] - T[i-1,j]**2 / (4.0*T[i-1,i-1])

        # Others
        for r in xrange(i,rr):
            for c in xrange(r+1,cc):
                Ti[r,c] = T[r,c] - T[i-1,r]*T[i-1,c] / (2*T[i-1,i-1])

        return Ti


    def prepare(self, rows, cols):
        r"""
        """
        # Unpack quadrature rules
        nodes = self._QR.get_nodes()
        weights = self._QR.get_weights()
        self._sqrtnodes = product(sqrt(nodes), axis=0)
        self._allweights = product(weights, axis=0)

        # Signs
        D = self._packet.get_dimension()
        x,y = indices((D, 2**D))
        self._allsigns = -2.0*(flipud((y>>x) & 1)-0.5)


    def do_nsd(self, row, col):
        r"""Evaluate a single integral :math:`\langle \phi_k | f(x) | \phi_l\rangle`
        by numerical steepest descent.

        :param row: The index :math:`k` of the basis function :math:`\phi_k` of :math:`\Phi`.
        :param col: The index :math:`l` of the basis function :math:`\phi^\prime_l` of :math:`\Phi^\prime`.
        :return: A single complex floating point number.
        """
        # TODO: Inline this function
        D = self._packet.get_dimension()

        T = self._oscillator["T"]
        U = self._oscillator["U"]
        X = self._oscillator["X"]
        b = self._oscillator["b"]

        # Compute all contour integrals
        results = []

        for i in xrange(2**D):
            # Signs
            signs = self._allsigns[:,i].reshape((D,))
            # Construct full paths
            ht = zeros_like(self._pathsqrtpart)
            for k in reversed(xrange(D)):
                ht[k,:] = ht[k,:] + signs[k] * self._pathsqrtpart[k,:]
                for j in xrange(k+1, D):
                    ht[k,:] = ht[k,:] - 0.5*T[k,j]/T[k,k] * ht[j,:]
            # Back transformation
            h = dot(conjugate(transpose(U)), ht) - dot(X, b)
            # Non-oscillatory parts
            fpath = self._envelope(h)
            opath = self._operator(h, entry=(row,col))
            #QQ = osign * pf * sum(fpath * dhts * sqrtnodes * allweights) / w**D
            QQ = sum(opath * fpath * self._quadrand)
            results.append(QQ)

        return self._prefactor * sum(squeeze(array(results)))


    def perform_quadrature(self, row, col):
        r"""Evaluates by numerical steepest descent the integral
        :math:`\langle \Phi_i | f | \Phi^\prime_j \rangle` for a polynomial
        function :math:`f(x)` with :math:`x \in \mathbb{R}^D`.

        :param row: The index :math:`i` of the component :math:`\Phi_i` of :math:`\Psi`.
        :param row: The index :math:`j` of the component :math:`\Phi^\prime_j` of :math:`\Psi^\prime`.
        :return: A single complex floating point number.
        """
        D = self._packet.get_dimension()

        # Combine oscillators
        Pibra = self._pacbra.get_parameters(key=("q", "p", "Q", "P"), component=row)
        Piket = self._packet.get_parameters(key=("q", "p", "Q", "P"), component=col)

        A, b, c = self.build_bilinear(Pibra, Piket)
        self._oscillator["A"] = A
        self._oscillator["b"] = b

        # Schur decomposition of A:
        #   A = U^H T U
        T, U = schur(A, output="complex")
        U = conjugate(transpose(U))
        self._oscillator["U"] = U

        # Oscillator updates
        for i in xrange(1, D):
            T = self.update_oscillator(T, i)
        self._oscillator["T"] = T

        # Take out diagonals
        Dk = diag(T).reshape((D,1))

        #
        X = inv(A + transpose(A))
        self._oscillator["X"] = X
        ctilde = c - 0.5 * dot(transpose(b), dot(X, b))

        # Prefactor
        eps = self._packet.get_eps()
        w = 1.0 / eps**2
        self._prefactor = exp(1.0j * w * ctilde)

        # Tau
        tk = self._nodes / w
        # Paths
        self._pathsqrtpart = sqrt(1.0j * tk / Dk)
        # Path derivatives
        dht = sqrt(1.0j / (4.0 * Dk * tk))
        dhtp = product(dht, axis=0)

        self._quadrand = dhtp * self._sqrtnodes * self._allweights / w**D


        # TODO: for all pairs <phi_k | ... | phi_l> do
        result = []
        # self._envelope = ...
        I = self.perform_nsd(row, col)
        result.append(I)

        return sum(I)


    def perform_build_matrix(self, row, col):
        raise NotImplementedError("")
