"""The WaveBlocks Project

This function makes a three dimensional stem plot.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from matplotlib.pyplot import gcf
import mpl_toolkits.mplot3d.art3d as art3d


def stem3d(u, v, w, fig=None, markerp="o"):
    r"""This function makes a three dimensional stem plot.
    """
    # Plot to the given axis instance or retrieve the current one
    if fig is None:
        fig = gcf()

    axes = fig.add_subplot(1, 1, 1, projection='3d')

    for ui, vi, wi in zip(u, v, w):
        line = art3d.Line3D(*zip((ui, vi, 0), (ui, vi, wi)), marker=markerp, markevery=(1, 1))
        axes.add_line(line)

    axes.set_xlim3d(u.min(), u.max())
    axes.set_ylim3d(v.min(), v.max())
    axes.set_zlim3d(w.min(), w.max())
