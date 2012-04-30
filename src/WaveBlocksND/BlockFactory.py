"""The WaveBlocks Project

This file contains a the block factory.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sympy

__all__ = ["BlockFactory"]


class BlockFactory(object):
    """A factory for :py:class:`Wavepacket` instances.
    """

    def __init__(self):
        pass


    # TODO: Consider "local" vs "global" description dicts
    # TODO: Consider putting defaults into "GlobalDefaults"


    def create_basis_shape(self, description):
        try:
            bs_type = description["type"]
        except:
            # Default setting
            bs_type = "HyperCubicShape"

        if bs_type == "HyperCubicShape":
            from HyperCubicShape import HyperCubicShape
            limits = description["limits"]
            BS = HyperCubicShape(limits)

        else:
            raise ValueError("Unknown basis shape type "+str(bs_type))

        return BS


    def create_wavepacket(self, description):

        wp_type = description["type"]

        if wp_type == "HagedornWavepacket":
            from HagedornWavepacket import HagedornWavepacket

            # Initialize a packet
            WP = HagedornWavepacket(description)

            # Set parameters
            if description.has_key("Pi"):
                Pi = description["Pi"]
                WP.set_parameters(Pi)

            # Configure basis shapes
            if description.has_key("basis_shapes"):
                for component, shapedescr in enumerate(description["basis_shapes"]):
                    BS = self.create_basis_shape(shapedescr)
                    WP.set_basis_shape(BS, component=component)

            # Set coefficients
            if description.has_key("coefficients"):
                for component, data in enumerate(description["coefficients"]):
                    for index, value in data:
                        WP.set_coefficient(component, index, value)

            # And the quadrature
            if description.has_key("quadrature"):
                QE = self.create_quadrature(description["quadrature"])
                WP.set_quadrature(QE)
            else:
                print("Warning: no quadrature specified!")

        else:
            raise ValueError("Unknown wavepacket type "+str(wb_type))

        return WP


    def create_quadrature(self, description):
        try:
            qe_type = description["type"]
        except:
            # Default setting
            qe_type = "InhomogeneousQuadrature"

        # TODO: Maybe denest QR initialization?
        if qe_type == "HomogeneousQuadrature":
            from HomogeneousQuadrature import HomogeneousQuadrature
            from GaussHermiteQR import GaussHermiteQR

            dim = description["dimension"]

            if dim == 1:
                order = description["qr_order"]
                assert type(order) == int
                QR = GaussHermiteQR(order)

            else:
                orders = description["qr_order"]
                # 'orders' should be a list here
                qrs_done = {}
                qrs = []
                for index, order in enumerate(orders):
                    if not qrs_done.has_key(order):
                        qrs_done[order] = index
                        qrs.append( GaussHermiteQR(order) )
                    else:
                        qrs.append( qrs[qrs_done[order]] )

                from TensorProductQR import TensorProductQR
                QR = TensorProductQR(qrs)

            QE = HomogeneousQuadrature(QR)

        elif qe_type == "InhomogeneousQuadrature":
            from InhomogeneousQuadrature import InhomogeneousQuadrature
            from GaussHermiteQR import GaussHermiteQR

            dim = description["dimension"]

            if dim == 1:
                order = description["qr_order"]
                assert type(order) == int
                QR = GaussHermiteQR(order)

            else:
                orders = description["qr_order"]
                # 'orders' should be a list here
                qrs_done = {}
                qrs = []
                for index, order in enumerate(orders):
                    if not qrs_done.has_key(order):
                        qrs_done[order] = index
                        qrs.append( GaussHermiteQR(order) )
                    else:
                        qrs.append( qrs[qrs_done[order]] )

                from TensorProductQR import TensorProductQR
                QR = TensorProductQR(qrs)

            QE = InhomogeneousQuadrature(QR)

        else:
            raise ValueError("Unknown basis shape type "+str(bs_type))

        return QE


    def create_quadrature_rule(self, description):
        qr_type = description["type"]

        if qr_type == "GaussHermiteQR":
            from GaussHermiteQR import GaussHermiteQR
            order = description["qr_order"]
            assert type(order) == int
            QR = GaussHermiteQR(order)

        elif qr_type == "TensorProductQR":
            from GaussHermiteQR import GaussHermiteQR
            orders = description["qr_order"]
            # 'orders' should be a list here
            qrs_done = {}
            qrs = []
            for index, order in enumerate(orders):
                if not qrs_done.has_key(order):
                    qrs_done[order] = index
                    qrs.append( GaussHermiteQR(order) )
                else:
                    qrs.append( qrs[qrs_done[order]] )

            from TensorProductQR import TensorProductQR
            QR = TensorProductQR(qrs)

        return QR
