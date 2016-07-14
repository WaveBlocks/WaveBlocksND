"""The WaveBlocks Project

This file contains the abstract base class for the computation of
the action of the gradient operator applied to a wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""

__all__ = ["Gradient"]


class Gradient(object):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x`.
    """

    def __init__(self):
        r"""
        :raise: :py:class:`NotImplementedError` This is an abstract base class.
        """
        raise NotImplementedError("'Gradient' is an abstract base class.")


    def apply_gradient(self, wavepacket, *, component=None):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x` on
        the basis functions :math:`\phi(x)` of a component :math:`\Phi_i` of a wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer or ``None``.
        :return: Extended basis shape :math:`\mathfrak{\dot{K}}` and new coefficients :math:`c^\prime`.
        """
        pass
