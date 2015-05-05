"""The WaveBlocks Project

Some selected functions for complex math.

@author: R. Bourquin
@copyright: Copyright (C) 2011, 2012, 2015 R. Bourquin
@license: Modified BSD License
"""

from numpy import array, hstack, cumsum, diff, around, abs, angle, exp, sqrt, pi, squeeze

__all__ = ["continuate", "cont_angle", "cont_sqrt"]


def continuate(data, jump=2.0*pi, reference=0.0):
    r"""Make the given data continuous by removing all jumps of size :math:`k\cdot\text{jump}`
    but not touching jumps of any other size. This can be used to overcome issues with
    the branch cut along the negative axis. There may be arise problems with jumps that
    are of size nearly :math:`k\cdot\text{jump}`.

    :param data: An array with the input data.
    :param jump: The basic size of jumps which will be removed. Default is :math:`2 \pi`.
    :param reference: This value allows the specify the starting point for continuation
                       explicitly. It can be used together with ``data``.
    :type reference: A single float number.
    """
    assert squeeze(reference).ndim == 0, "The 'reference' must be a scalar value."
    data = hstack([array(reference), array(data).reshape(-1)])
    offsets = cumsum(around(diff(data)/(1.0*jump)))
    return (data - jump*hstack([0.0, offsets]))[1:]


def cont_angle(data, reference=None):
    r"""Compute the angle of a complex number *not* constrained to
    the principal value and avoiding discontinuities at the branch cut.
    This function just applies 'continuate(.)' to the complex phase.

    :param data: An array with the input data.
    :param reference: This value allows the specify the starting point for
                      continuation explicitly. It can be used together with ``data``.
    :type reference: A single float number.
    """
    if reference is None:
        # Return just cont_f(x)
        return continuate(angle(data))
    else:
        # Return a 2-tuple (cont_f(x), new_reference)
        result = continuate(angle(data), reference=reference)
        reference = result
        return result, reference


def cont_sqrt(data, reference=None):
    r"""Compute the complex square root (following the Riemann surface)
    yields a result *not* constrained to the principal value and avoiding
    discontinuities at the branch cut. This function applies 'continuate(.)'
    to the complex phase and computes for
    :math:`z = r \exp \left(i \phi \right)` the complex square root
    its square root according to the formula
    :math:`\sqrt{z} = \sqrt{r} \exp \left(i \frac{\phi}{2} \right)`.

    :param data: An array with the input data.
    :param reference: This value allows the specify the starting point for
                      continuation explicitly. It can be used together with ``data``.
    :type reference: A single float number.
    """
    if reference is None:
        # Return just cont_f(x)
        return sqrt(abs(data)) * exp(0.5j*continuate(angle(data)))
    else:
        # Return a 2-tuple (cont_f(x), new_reference)
        phi = continuate(angle(data), reference=reference)
        result = sqrt(abs(data)) * exp(0.5j*phi)
        reference = phi[0]
        # TODO: Rethink what 'reference' to return for an array of input values.
        return result, reference


class ContinuousSqrt(object):
    r"""Class for computing continuous square roots.
    All implementation details about referencing is
    hidden from the user. The class is side-effect free.
    """

    def __init__(self, reference=0.0j):
        self.set(reference)

    def clone(self):
        return ContinuousSqrt(reference=self._reference)

    def __call__(self, radicand):
        # Updating the reference is idempotent for identical radicands.
        root, reference = cont_sqrt(radicand, reference=self._reference)
        self._reference = reference
        return root

    def set(self, reference):
        assert squeeze(reference).ndim == 0, "The 'reference' must be a scalar value."
        self._reference = array(reference)

    def get(self):
        return self._reference.copy()
