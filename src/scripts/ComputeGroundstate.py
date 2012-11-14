"""The WaveBlocks Project

This file contains the class which represents a homogeneous Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from numpy import *
from scipy import *
from numpy import squeeze, real, zeros, zeros_like, conj, sum
from scipy.optimize import fmin
from scipy.linalg import sqrtm, inv, eigh

from WaveBlocksND import BlockFactory
from WaveBlocksND import ParameterLoader
from WaveBlocksND import ParameterProvider
from WaveBlocksND import GradientHAWP
from WaveBlocksND import IOManager


def compute_groundstate(params):
    r"""
    """
    D = params["dimension"]
    N = params["groundstate_of_level"]

    # Create output file now, in case this fails we did not waste computations
    IOM = IOManager()
    IOM.create_file(params, filename="groundstate.hdf5")
    gid = IOM.create_group()
    IOM.create_block(groupid=gid)

    BF = BlockFactory()
    # Create the potential
    V = BF.create_potential(params)
    V.calculate_local_quadratic()

    # Minimize the potential to find q0
    f = lambda x: real((squeeze(V.evaluate_at(x)[N])))
    # Start with an offset because exact 0.0 values can give
    # issues, especially with the Hessian evaluation. This way
    # the minimizer will always stay away from zero a tiny bit.
    x0 = 2.0*ones(D)
    q0 = fmin(f, x0, xtol=1e-12)
    q0 = q0.reshape((D,1))

    # We are at the minimum with no momentum
    p0 = zeros_like(q0)

    # Compute spreads now
    # Q_0 = H^(-1/4)
    H = V.evaluate_hessian_at(q0)
    Q0 = inv(sqrtm(sqrtm(H)))
    # Take P_00 = i Q_0^(-1)
    P0 = 1.0j * inv(Q0)

    #
    print("Parameter values are:")
    print(" q0:")
    print(str(q0))
    print(" p0:")
    print(str(p0))
    print(" Q0:")
    print(str(Q0))
    print(" P0:")
    print(str(P0))
    # Consistency check
    print(" consistency:")
    print(str(conj(Q0)*P0 - conj(P0)*Q0))

    # Next find the new coefficients c'
    HAWP = BF.create_wavepacket(params["hawp_template"])

    # Set the parameter values
    Pi = HAWP.get_parameters()
    Pi[0] = q0
    Pi[1] = p0
    Pi[2] = Q0
    Pi[3] = P0
    HAWP.set_parameters(Pi)

    # Next compute the matrix M_ij = <phi_i | T + V | phi_j>
    # First the potential part
    HQ = BF.create_quadrature(params["quadrature"])

    opV = lambda x, q: V.evaluate_at(x)
    MV = HQ.build_matrix(HAWP, operator=opV)

    # The kinetic part
    MT = zeros_like(MV, dtype=complexfloating)
    GR = GradientHAWP()
    BS = HAWP.get_basis_shapes(N)

    vects = {}
    for i in BS:
        z = zeros_like(HAWP.get_coefficient_vector(), dtype=complexfloating)
        HAWP.set_coefficient_vector(z)
        HAWP.set_coefficient(N, i, 1.0)
        Kn, cnew = GR.apply_gradient(HAWP, N)
        vects[i] = cnew

    for j in BS:
        for k in BS:
            cj = vects[j]
            ck = vects[k]
            entry = 0.5 * squeeze(sum(conj(cj) * ck))
            MT[BS[j], BS[k]] = entry

    # Find eigenvalues and eigenvectors of the whole matrix
    M = MT + MV
    ew, ev = eigh(M)
    ind = argsort(ew)
    emin_ind = ind[0]
    cfinal = ev[:, emin_ind]
    ew[emin_ind]

    cfinal = cfinal.reshape((-1, 1))
    HAWP.set_coefficient_vector(cfinal)

    # Save all the wavepacket data
    IOM.add_wavepacket(params)
    IOM.save_wavepacket_description(HAWP.get_description())
    for shape in HAWP.get_basis_shapes():
        IOM.save_wavepacket_basisshapes(shape)
    IOM.save_wavepacket_parameters(HAWP.get_parameters(), timestep=0)
    IOM.save_wavepacket_coefficients(HAWP.get_coefficients(), HAWP.get_basis_shapes(), timestep=0)

    IOM.finalize()




if __name__ == "__main__":


    # Read the path for the configuration file we use for this simulation.
    try:
        parametersfile = sys.argv[1]
    except IndexError:
        raise ValueError("No configuration file given!")

    print("Using configuration from file: " + parametersfile)

    # Set up the parameter provider singleton
    PA = ParameterLoader().load_from_file(parametersfile)
    compute_groundstate(PA)

    print("Groundstate computation finished")
