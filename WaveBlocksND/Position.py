"""The WaveBlocks Project

This file contains the abstract base class for the computation of
the action of the position operator applied to a wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

__all__ = ["Position"]


class Position(object):
    r"""This class implements the computation of the action of the
    position operator :math:`x`.
    """

    def __init__(self):
        r"""
        :raise: :py:class:`NotImplementedError` This is an abstract base class.
        """
        raise NotImplementedError("'Position' is an abstract base class.")


    def apply_position(self, wavepacket, *, component=None):
        r"""Compute the effect of the position operator :math:`x` on the basis
        functions :math:`\phi(x)` of a component :math:`\Phi_i` of a wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer or ``None``.
        :return: Extended basis shape :math:`\mathfrak{\dot{K}}` and new coefficients :math:`c^\prime`.
        """
        pass
