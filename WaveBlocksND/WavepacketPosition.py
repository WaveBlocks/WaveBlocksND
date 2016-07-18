"""The WaveBlocks Project

Compute the action of the position operator applied to a wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND.Position import Position

__all__ = ["WavepacketPosition"]


class WavepacketPosition(Position):
    r"""This class implements the computation of the action of the position
    operator :math:`x` applied to a Hagedorn wavepacket :math:`\Psi`.
    """

    def __init__(self):
        r"""
        """
        pass


    def apply_position(self, wavepacket, *, component=None, as_packet=True):
        r"""Compute the effect of the position operator :math:`x` on the basis
        functions :math:`\phi(x)` of a component :math:`\Phi_i` of the Hagedorn
        wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer or ``None``.
        :param as_packet: Whether to return a full packet.
        :type as_packet: Boolean, default is ``True``.
        :return: A list of length :math:`N` or a single pair of extended basis shapes :math:`\mathfrak{\dot{K}}_i`
                 and new coefficients :math:`c^\prime_i`.
                 If requested, copies of the original wavepacket are returned with these new values set.
                 There are :math:`D` packets, one for each space variable component :math:`x_d`
                 of the position operator.
        """
        if component is not None:
            components = [component]
        else:
            N = wavepacket.get_number_components()
            components = range(N)

        # Apply the position operator to all components
        positionparts = [self.apply_position_component(wavepacket, n) for n in components]

        if as_packet is True:
            D = wavepacket.get_dimension()
            new_wps = []
            for d in range(D):
                hawp_new = wavepacket.clone()
                for i, component in enumerate(components):
                    Ke, cnew = positionparts[i]
                    hawp_new.set_basis_shapes(Ke, component=component)
                    hawp_new.set_coefficients(cnew[:, d], component=component)
                new_wps.append(hawp_new)
            return new_wps
        else:
            if component is not None:
                return positionparts[0]
            return positionparts
