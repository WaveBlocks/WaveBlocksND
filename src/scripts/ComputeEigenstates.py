"""The WaveBlocks Project

This file contains code for computing the eigenstates of a given
potential in terms of Hagedorn wavepackets.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014 R. Bourquin
@license: Modified BSD License
"""

import argparse
from numpy import (argsort, atleast_1d, atleast_2d, complexfloating, conjugate, dot, ones, zeros,
                   real, identity, squeeze, sum, transpose, zeros_like, argmax, angle, abs, pi)
from scipy.optimize import fmin
from scipy.linalg import sqrtm, inv, eigh, norm

from WaveBlocksND import BlockFactory
from WaveBlocksND import GradientHAWP
from WaveBlocksND import IOManager
from WaveBlocksND import ParameterLoader


def compute_eigenstate(parameters, filename="eigenstates.hdf5", computepq=True, computePQ=True):
    r"""
    Special variables necessary in configuration:

    * eigenstate_of_level (default: 0)
    * eigenstates_indices (default: [0])
    * starting_point (default: (2, ..., 2))
    * hawp_template
    * innerproduct
    """
    D = parameters["dimension"]

    if "eigenstate_of_level" in parameters:
        N = parameters["eigenstate_of_level"]
    else:
        # Upper-most potential surface
        N = 0

    # Create output file now, in case this fails we did not waste computation time
    IOM = IOManager()
    IOM.create_file(filename)

    # Save the simulation parameters
    IOM.add_parameters()
    IOM.save_parameters(parameters)

    gid = IOM.create_group()

    BF = BlockFactory()
    # Create the potential
    V = BF.create_potential(parameters)
    V.calculate_local_quadratic()

    # Compute position and momentum
    if computepq:
        # Minimize the potential to find q0
        f = lambda x: real((squeeze(V.evaluate_at(x)[N])))
        # Start with an offset because exact 0.0 values can give
        # issues, especially with the Hessian evaluation. This way
        # the minimizer will always stay away from zero a tiny bit.
        # The current starting point can give issues if the potential
        # is stationary at the point (2, ..., 2) but that is less likely.
        if "starting_point" in parameters:
            x0 = atleast_1d(parameters["starting_point"])
        else:
            x0 = 0.5 * ones(D)

        q0 = fmin(f, x0, xtol=1e-12)
        q0 = q0.reshape((D,1))

        # We are at the minimum with no momentum
        p0 = zeros_like(q0)
    else:
        if "q0" in parameters:
            q0 = atleast_2d(parameters["q0"])
        else:
            q0 = zeros((D,1))
        if "p0" in parameters:
            p0 = atleast_2d(parameters["p0"])
        else:
            p0 = zeros((D,1))

    # Compute spreads
    if computePQ:
        # Q_0 = H^(-1/4)
        H = V.evaluate_hessian_at(q0)
        Q0 = inv(sqrtm(sqrtm(H)))
        # P_0 = i Q_0^(-1)
        P0 = 1.0j * inv(Q0)
    else:
        if "Q0" in parameters:
            Q0 = atleast_2d(parameters["Q0"])
        else:
            Q0 = identity(D)
        if "P0" in parameters:
            P0 = atleast_2d(parameters["P0"])
        else:
            P0 = 1.0j * inv(Q0)

    # The parameter set Pi
    print(70*"-")
    print("Parameter values are:")
    print("---------------------")
    print(" q0:")
    print(str(q0))
    print(" p0:")
    print(str(p0))
    print(" Q0:")
    print(str(Q0))
    print(" P0:")
    print(str(P0))
    # Consistency check
    print(" Consistency check:")
    print("   P^T Q - Q^T P  =?=  0")
    print(dot(P0.T,Q0) - dot(Q0.T,P0))
    print("   Q^H P - P^H Q  =?=  2i")
    print(dot(transpose(conjugate(Q0)), P0) - dot(transpose(conjugate(P0)), Q0))

    # Next find the new coefficients c'
    HAWP = BF.create_wavepacket(parameters["hawp_template"])

    # Set the parameter values
    Pi = HAWP.get_parameters()
    Pi[0] = q0
    Pi[1] = p0
    Pi[2] = Q0
    Pi[3] = P0
    HAWP.set_parameters(Pi)

    # Next compute the matrix M_ij = <phi_i | T + V | phi_j>
    # The potential part
    HQ = BF.create_inner_product(parameters["innerproduct"])

    opV = lambda x, q, entry: V.evaluate_at(x, entry=entry)
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
        Kn, cnew = GR.apply_gradient(HAWP, component=N, as_packet=False)
        vects[i] = cnew

    for j in BS:
        for k in BS:
            cj = vects[j]
            ck = vects[k]
            entry = 0.5 * squeeze(sum(conjugate(cj) * ck))
            MT[BS[j], BS[k]] = entry

    # Find eigenvalues and eigenvectors of the whole matrix
    M = MT + MV
    ew, ev = eigh(M)
    ind = argsort(ew)

    # Build the requested energy levels and states
    if "eigenstates_indices" in parameters:
        states = parameters["eigenstates_indices"]
    else:
        # Groundstate only
        states = [0]

    BS = HAWP.get_basis_shapes(component=0)

    KEY = ("q","p","Q","P","S","adQ")

    print(70*"-")
    for state in states:
        if state > BS.get_basis_size():
            print("Warning: can not compute energy level "+state+" with basis size of "+str(BS))
            continue

        index = ind[state]

        coeffs = ev[:,index]
        energy = ew[index]

        # Try to resolve ambiguities in sign
        imax = argmax(abs(coeffs))
        a = abs(angle(coeffs[imax]))
        if a > pi/2.0:
            coeffs *= -1

        print("State: "+str(state))
        print("Energy: "+str(energy))
        print("Coefficients: \n")
        print(str(coeffs))
        print(70*"-")

        HAWP.set_coefficient_vector(coeffs.reshape((-1, 1)))

        # Save all the wavepacket data
        bid = IOM.create_block(groupid=gid)
        IOM.add_wavepacket(parameters, blockid=bid, key=KEY)
        IOM.save_wavepacket(HAWP, 0, blockid=bid, key=KEY)

    IOM.finalize()

    if norm(q0) > 1000:
        print("+----------------------------------+")
        print("| Run-away minimum?                |")
        print("| Maybe try different:             |")
        print("|   starting_point = [x0, y0, ...] |")
        print("+----------------------------------+")




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("parametersfile",
                        type = str,
                        help = "The configuration parameters file")

    parser.add_argument("-o", "--outputfile",
                        type = str,
                        help = "The output data file",
                        nargs = "?",
                        default = "eigenstates.hdf5")

    parser.add_argument("--nopq",
                        help = "Do not compute the parameters q and p",
                        action = "store_false")

    parser.add_argument("--noPQ",
                        help = "Do not compute the parameters Q and P",
                        action = "store_false")

    args = parser.parse_args()

    # Read the path for the configuration file we use for this simulation.
    print("Using configuration from file: " + args.parametersfile)

    # Set up the parameter provider singleton
    PA = ParameterLoader().load_from_file(args.parametersfile)
    compute_eigenstate(PA, filename=args.outputfile, computepq=args.nopq, computePQ=args.noPQ)

    print("Eigenstate computation finished")
