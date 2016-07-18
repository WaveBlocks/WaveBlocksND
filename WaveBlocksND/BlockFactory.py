
"""The WaveBlocks Project

This file contains a the block factory.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

from copy import deepcopy

__all__ = ["BlockFactory"]


class BlockFactory(object):
    """A factory to create instances of various classes
    based on a simple description ``dict``.
    """

    def __init__(self):

        # Load different factory methods
        from WaveBlocksND import GridFactory
        self.__dict__["create_grid"] = GridFactory.create_grid

        from WaveBlocksND import PotentialFactory
        self.__dict__["create_potential"] = PotentialFactory.create_potential

        from WaveBlocksND import MatrixExponentialFactory
        self.__dict__["create_matrixexponential"] = MatrixExponentialFactory.create_matrixexponential


    # TODO: Consider "local" vs "global" description dicts
    # TODO: Consider putting defaults into "GlobalDefaults"


    def create_basis_shape(self, description):
        try:
            bs_type = description["type"]
        except:
            # Default setting
            bs_type = "HyperCubicShape"

        if bs_type == "HyperCubicShape":
            from WaveBlocksND.HyperCubicShape import HyperCubicShape
            limits = description["limits"]
            BS = HyperCubicShape(limits)

        elif bs_type == "SimplexShape":
            from WaveBlocksND.SimplexShape import SimplexShape
            K = description["K"]
            D = description["dimension"]
            BS = SimplexShape(D, K)

        elif bs_type == "HyperbolicCutShape":
            from WaveBlocksND.HyperbolicCutShape import HyperbolicCutShape
            K = description["K"]
            D = description["dimension"]
            BS = HyperbolicCutShape(D, K)

        else:
            raise ValueError("Unknown basis shape type {}".format(bs_type))

        return BS


    def create_wavepacket(self, description):

        wp_type = description["type"]

        if wp_type == "HagedornWavepacket":
            from WaveBlocksND.HagedornWavepacket import HagedornWavepacket

            # Initialize a packet
            WP = HagedornWavepacket(description["dimension"],
                                    description["ncomponents"],
                                    description["eps"])

            # Set parameters
            if "Pi" in description:
                Pi = description["Pi"]
                WP.set_parameters(Pi)

            # Configure basis shapes
            if "basis_shapes" in description:
                for component, shapedescr in enumerate(description["basis_shapes"]):
                    BS = self.create_basis_shape(shapedescr)
                    WP.set_basis_shapes(BS, component=component)

            # Set coefficients
            if "coefficients" in description:
                for component, data in enumerate(description["coefficients"]):
                    BS = WP.get_basis_shapes(component=component)
                    for index, value in data:
                        if BS.contains(index):
                            WP.set_coefficient(component, index, value)
                        else:
                            print("Warning: dropped coefficient with index {}".format(index))

            # And the inner product
            if "innerproduct" in description:
                IP = self.create_inner_product(description["innerproduct"])
                WP.set_innerproduct(IP)
            else:
                print("Warning: no inner product specified!")

        elif wp_type == "HagedornWavepacketInhomogeneous":
            from WaveBlocksND.HagedornWavepacketInhomogeneous import HagedornWavepacketInhomogeneous

            # Initialize a packet
            WP = HagedornWavepacketInhomogeneous(description["dimension"],
                                                 description["ncomponents"],
                                                 description["eps"])

            # Set parameters
            if "Pi" in description:
                Pi = description["Pi"]
                WP.set_parameters(Pi)

            # Configure basis shapes
            if "basis_shapes" in description:
                for component, shapedescr in enumerate(description["basis_shapes"]):
                    BS = self.create_basis_shape(shapedescr)
                    WP.set_basis_shapes(BS, component=component)

            # Set coefficients
            if "coefficients" in description:
                for component, data in enumerate(description["coefficients"]):
                    for index, value in data:
                        WP.set_coefficient(component, index, value)

            # And the quadrature
            if "innerproduct" in description:
                IP = self.create_inner_product(description["innerproduct"])
                WP.set_innerproduct(IP)
            else:
                print("Warning: no inner product specified!")

        else:
            raise ValueError("Unknown wavepacket type {}".format(wp_type))

        return WP


    def create_inner_product(self, description):
        #
        try:
            ip_type = description["type"]
        except:
            # Default setting
            print("Warning: Using fall-back inner product of type 'InhomogeneousInnerProduct'!")
            # TODO: Maybe change to homogeneous one?
            ip_type = "InhomogeneousInnerProduct"

        if ip_type == "HomogeneousInnerProduct":
            from WaveBlocksND.HomogeneousInnerProduct import HomogeneousInnerProduct
            QE = self.create_quadrature(description["delegate"])
            IP = HomogeneousInnerProduct(QE)

        elif ip_type == "InhomogeneousInnerProduct":
            from WaveBlocksND.InhomogeneousInnerProduct import InhomogeneousInnerProduct
            QE = self.create_quadrature(description["delegate"])
            IP = InhomogeneousInnerProduct(QE)

        else:
            raise ValueError("Unknown inner product type {}".format(ip_type))

        return IP


    def create_quadrature(self, description):
        # NOTE: The difference between Quadrature and InnerProduct here
        #       is for backward compatibility but will be removed soon.
        try:
            qe_type = description["type"]
        except:
            # Default setting
            qe_type = "DirectInhomogeneousQuadrature"

        # TODO: Maybe denest QR initialization?
        if qe_type == "DirectHomogeneousQuadrature":
            from WaveBlocksND.DirectHomogeneousQuadrature import DirectHomogeneousQuadrature
            QR = self.create_quadrature_rule(description["qr"])
            QE = DirectHomogeneousQuadrature(QR)

        elif qe_type == "DirectInhomogeneousQuadrature":
            from WaveBlocksND.DirectInhomogeneousQuadrature import DirectInhomogeneousQuadrature
            QR = self.create_quadrature_rule(description["qr"])
            QE = DirectInhomogeneousQuadrature(QR)

        elif qe_type == "NSDInhomogeneous":
            from WaveBlocksND.NSDInhomogeneous import NSDInhomogeneous
            QR = self.create_quadrature_rule(description["qr"])
            QE = NSDInhomogeneous(QR)

        elif qe_type == "SymbolicIntegral":
            from WaveBlocksND.SymbolicIntegral import SymbolicIntegral
            QE = SymbolicIntegral()

        else:
            raise ValueError("Unknown quadrature type {}".format(qe_type))

        return QE


    def create_quadrature_rule(self, description):
        qr_type = description["type"]

        if "options" in description:
            op = deepcopy(description["options"])
        else:
            # Per default, adapt qr to follow dynamics
            op = {"transform": True}

        if qr_type == "GaussHermiteQR":
            from WaveBlocksND.GaussHermiteQR import GaussHermiteQR
            order = description["order"]
            assert type(order) == int
            QR = GaussHermiteQR(order, options=op)

        elif qr_type == "GaussHermiteOriginalQR":
            from WaveBlocksND.GaussHermiteOriginalQR import GaussHermiteOriginalQR
            order = description["order"]
            assert type(order) == int
            QR = GaussHermiteOriginalQR(order, options=op)

        elif qr_type == "TrapezoidalQR":
            from WaveBlocksND.TrapezoidalQR import TrapezoidalQR
            left = description["left"]
            right = description["right"]
            order = description["order"]
            assert type(order) == int
            QR = TrapezoidalQR(left, right, order, options=op)

        elif qr_type == "GaussLaguerreQR":
            from WaveBlocksND.GaussLaguerreQR import GaussLaguerreQR
            order = description["order"]
            a = description["a"]
            assert type(order) == int
            QR = GaussLaguerreQR(order, a=a, options=op)

        elif qr_type == "TensorProductQR":
            from WaveBlocksND.TensorProductQR import TensorProductQR
            # Iteratively create all quadrature rules necessary
            qrs = [self.create_quadrature_rule(desc) for desc in description["qr_rules"]]
            QR = TensorProductQR(qrs, options=op)

        else:
            raise ValueError("Unknown quadrature rule type {}".format(qr_type))

        return QR


    def create_propagator(self, description, *args, **kwargs):
        prop_type = description["propagator"]

        if prop_type == "fourier":
            from WaveBlocksND.FourierPropagator import FourierPropagator
            propagator = FourierPropagator(description, *args, **kwargs)

        elif prop_type == "hagedorn":
            from WaveBlocksND.HagedornPropagator import HagedornPropagator
            propagator = HagedornPropagator(description, *args, **kwargs)

        elif prop_type == "semiclassical":
            from WaveBlocksND.SemiclassicalPropagator import SemiclassicalPropagator
            propagator = SemiclassicalPropagator(description, *args, **kwargs)

        elif prop_type == "magnus_split":
            from WaveBlocksND.MagnusPropagator import MagnusPropagator
            propagator = MagnusPropagator(description, *args, **kwargs)

        elif prop_type == "McL42sc":
            from WaveBlocksND.McL42scPropagator import McL42scPropagator
            propagator = McL42scPropagator(description, *args, **kwargs)

        elif prop_type == "McL84sc":
            from WaveBlocksND.McL84scPropagator import McL84scPropagator
            propagator = McL84scPropagator(description, *args, **kwargs)

        elif prop_type == "Pre764sc":
            from WaveBlocksND.Pre764scPropagator import Pre764scPropagator
            propagator = Pre764scPropagator(description, *args, **kwargs)

        elif prop_type == "hagedorn_inhomog":
            from WaveBlocksND.HagedornPropagatorInhomogeneous import HagedornPropagatorInhomogeneous
            propagator = HagedornPropagatorInhomogeneous(description, *args, **kwargs)

        elif prop_type == "hagedorn_psi":
            from WaveBlocksND.HagedornPropagatorPsi import HagedornPropagatorPsi
            propagator = HagedornPropagatorPsi(description, *args, **kwargs)

        else:
            raise ValueError("Unknown propagator type {}".format(prop_type))

        return propagator
