"""The WaveBlocks Project

This file contains the block factory.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2013, 2014, 2015 R. Bourquin
@license: Modified BSD License
"""

import importlib
from functools import partial
from copy import deepcopy

__all__ = ["BlockFactory"]


class BlockFactory(object):
    """A factory to create instances of various classes
    based on a simple description ``dict``.
    """

    def __init__(self):

        # Load different factory methods
        import PotentialFactory
        self.__dict__["create_potential"] = PotentialFactory.create_potential

    # TODO: Consider "local" vs "global" description dicts
    # TODO: Consider putting defaults into "GlobalDefaults"


    def create_matrixexponential(self, description):
        """Returns the requested matrix exponential routine.

        :param description: A :py:class:`ParameterProvider` instance containing at least the
                            key ``matrix_exponential`` and depending on its values more keys.
        """
        method = description["matrix_exponential"]

        if method == "pade":
            from MatrixExponential import matrix_exp_pade
            return matrix_exp_pade
        elif method == "arnoldi":
            import GlobalDefaults
            from MatrixExponential import matrix_exp_arnoldi
            arnoldi_steps = description.get("arnoldi_steps", GlobalDefaults.arnoldi_steps)
            return partial(matrix_exp_arnoldi, k=arnoldi_steps)
        else:
            raise ValueError("Unknown matrix exponential algorithm")


    def create_grid(self, description):
        """The method that creates a :py:class:`Grid` instance and decides
        which subclass to instantiate depending on the given description.

        :param description: A ``description`` (``dict`` or :py:class:`ParameterProvider` instance)
                            with all necessary parameters.
        :return: An adequate :py:class:`Grid` instance.
        """
        grid_type = description.get("type", "TensorProductGrid")

        if grid_type == "TensorProductGrid":
            from . import TensorProductGrid
            limits = description["limits"]
            number_nodes = description["number_nodes"]

            # TODO: Improve for one|multiple values: limits = "D*(a,b)" || "[(a1,b1), (a2,b2), ...]"

            grid = TensorProductGrid(limits, number_nodes)

        return grid


    def create_basis_shape(self, description):
        bs_args = description.copy()
        bs_type = bs_args.pop("type", "HyperCubicShape")

        try:
            M = importlib.import_module(bs_type, package="WaveBlocksND")
        except ImportError:
            raise ValueError("Unknown basis shape type "+str(bs_type))

        C = getattr(M, bs_type)
        BS = C(**bs_args)

        return BS


    def create_wavepacket(self, description):

        wp_type = description["type"]

        if wp_type == "HagedornWavepacket":
            from . import HagedornWavepacket

            # Initialize a packet
            WP = HagedornWavepacket(description["dimension"],
                                    description["ncomponents"],
                                    description["eps"])

            # Set parameters
            if description.has_key("Pi"):
                Pi = description["Pi"]
                WP.set_parameters(Pi)

            # Configure basis shapes
            if description.has_key("basis_shapes"):
                for component, shapedescr in enumerate(description["basis_shapes"]):
                    BS = self.create_basis_shape(shapedescr)
                    WP.set_basis_shapes(BS, component=component)

            # Set coefficients
            if description.has_key("coefficients"):
                for component, data in enumerate(description["coefficients"]):
                    BS = WP.get_basis_shapes(component=component)
                    for index, value in data:
                        if BS.contains(index):
                            WP.set_coefficient(component, index, value)
                        else:
                            print("Warning: dropped coefficient with index "+str(index))

            # And the inner product
            if description.has_key("innerproduct"):
                IP = self.create_inner_product(description["innerproduct"])
                WP.set_innerproduct(IP)
            else:
                print("Warning: no inner product specified!")

        elif wp_type == "HagedornWavepacketInhomogeneous":
            from . import HagedornWavepacketInhomogeneous

            # Initialize a packet
            WP = HagedornWavepacketInhomogeneous(description["dimension"],
                                    description["ncomponents"],
                                    description["eps"])

            # Set parameters
            if description.has_key("Pi"):
                Pi = description["Pi"]
                WP.set_parameters(Pi)

            # Configure basis shapes
            if description.has_key("basis_shapes"):
                for component, shapedescr in enumerate(description["basis_shapes"]):
                    BS = self.create_basis_shape(shapedescr)
                    WP.set_basis_shapes(BS, component=component)

            # Set coefficients
            if description.has_key("coefficients"):
                for component, data in enumerate(description["coefficients"]):
                    for index, value in data:
                        WP.set_coefficient(component, index, value)

            # And the quadrature
            if description.has_key("innerproduct"):
                IP = self.create_inner_product(description["innerproduct"])
                WP.set_innerproduct(IP)
            else:
                print("Warning: no inner product specified!")

        else:
            raise ValueError("Unknown wavepacket type "+str(wp_type))

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
            from . import HomogeneousInnerProduct
            QE = self.create_quadrature(description["delegate"])
            IP = HomogeneousInnerProduct(QE)

        elif ip_type == "InhomogeneousInnerProduct":
            from . import InhomogeneousInnerProduct
            QE = self.create_quadrature(description["delegate"])
            IP = InhomogeneousInnerProduct(QE)

        else:
            raise ValueError("Unknown inner product type "+str(ip_type))

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
            from . import DirectHomogeneousQuadrature
            QR = self.create_quadrature_rule(description["qr"])
            QE = DirectHomogeneousQuadrature(QR)

        elif qe_type == "DirectInhomogeneousQuadrature":
            from . import DirectInhomogeneousQuadrature
            QR = self.create_quadrature_rule(description["qr"])
            QE = DirectInhomogeneousQuadrature(QR)

        elif qe_type == "NSDInhomogeneous":
            from . import NSDInhomogeneous
            QR = self.create_quadrature_rule(description["qr"])
            QE = NSDInhomogeneous(QR)

        elif qe_type == "SymbolicIntegral":
            from . import SymbolicIntegral
            QE = SymbolicIntegral()

        else:
            raise ValueError("Unknown quadrature type "+str(qe_type))

        return QE


    def create_quadrature_rule(self, description):
        qr_type = description["type"]

        if description.has_key("options"):
            op = deepcopy(description["options"])
        else:
            # Per default, adapt qr to follow dynamics
            op = {"transform":True}

        if qr_type == "GaussHermiteQR":
            from . import GaussHermiteQR
            order = description["order"]
            assert type(order) == int
            QR = GaussHermiteQR(order, options=op)

        elif qr_type == "GaussHermiteOriginalQR":
            from . import GaussHermiteOriginalQR
            order = description["order"]
            assert type(order) == int
            QR = GaussHermiteOriginalQR(order, options=op)

        elif qr_type == "TrapezoidalQR":
            from . import TrapezoidalQR
            left = description["left"]
            right = description["right"]
            order = description["order"]
            assert type(order) == int
            QR = TrapezoidalQR(left, right, order, options=op)

        elif qr_type == "GaussLaguerreQR":
            from . import GaussLaguerreQR
            order = description["order"]
            a = description["a"]
            assert type(order) == int
            QR = GaussLaguerreQR(order, a=a, options=op)

        elif qr_type == "TensorProductQR":
            from . import TensorProductQR
            # Iteratively create all quadrature rules necessary
            qrs = [ self.create_quadrature_rule(desc) for desc in description["qr_rules"] ]
            QR = TensorProductQR(qrs, options=op)

        else:
            raise ValueError("Unknown quadrature rule type "+str(qr_type))

        return QR
