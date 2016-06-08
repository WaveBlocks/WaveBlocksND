"""The WaveBlocks Project

Compute the action of the gradient operator applied to a wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

from WaveBlocksND.Gradient import Gradient

__all__ = ["WavepacketGradient"]


class WavepacketGradient(Gradient):
    r"""This class implements the computation of the action of the
    gradient operator :math:`-i \varepsilon^2 \nabla_x` applied to
    a Hagedorn wavepacket :math:`\Psi`.
    """

    def __init__(self):
        r"""
        """
        pass


    def apply_gradient(self, wavepacket, *, component=None, as_packet=True):
        r"""Compute the effect of the gradient operator :math:`-i \varepsilon^2 \nabla_x` on the basis
        functions :math:`\phi(x)` of a component :math:`\Phi_i` of the Hagedorn wavepacket :math:`\Psi`.

        :param wavepacket: The wavepacket :math:`\Psi` containing :math:`\Phi_i`.
        :type wavepacket: A :py:class:`HagedornWavepacketBase` subclass instance.
        :param component: The index :math:`i` of the component :math:`\Phi_i`.
        :type component: Integer or ``None``.
        :param as_packet: Whether to return a full packet.
        :type as_packet: Boolean, default is ``True``.
        :return: A list of length :math:`N` or a single pair of extended basis shapes :math:`\mathfrak{\dot{K}}_i`
                 and new coefficients :math:`c^\prime_i`.
                 If requested, copies of the original wavepacket are returned with these new values set.
                 There are :math:`D` packets, one for each space variable component :math:`\partial_{x_d}`
                 of the gradient.
        """
        if component is not None:
            components = [component]
        else:
            N = wavepacket.get_number_components()
            components = range(N)

        # Compute the gradients of all components
        gradients = [self.apply_gradient_component(wavepacket, n) for n in components]

        if as_packet is True:
            D = wavepacket.get_dimension()
            new_wps = []
            for d in range(D):
                hawp_new = wavepacket.clone()
                for i, component in enumerate(components):
                    Ke, cnew = gradients[i]
                    hawp_new.set_basis_shapes(Ke, component=component)
                    hawp_new.set_coefficients(cnew[:, d], component=component)
                new_wps.append(hawp_new)
            return new_wps
        else:
            if component is not None:
                return gradients[0]
            return gradients
