"""The WaveBlocks Project

This file contains unit tests for the
TensorProductGrid class.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import array

from WaveBlocksND import TensorProductGrid


class TestTensorProductGrid:

    def build_grid(self):
        # Params
        d = {}
        d['dimension'] = 2
        d['grid_limits'] = [(-1, 2),(0,3)]
        d['grid_number_nodes'] = [8,10]

        Grid = TensorProductGrid(d)

        return Grid


    def test_dimension(self):
        G = self.build_grid()
        d = G.get_dimension()
        assert d == 2


    def test_limits(self):
        G = self.build_grid()

        l = G.get_limits()
        assert l == [(-1, 2), (0, 3)]

        l = G.get_limits(axes=0)
        assert l == [(-1, 2)]

        l = G.get_limits(axes=[1])
        assert l == [(0, 3)]

        l = G.get_limits(axes=[0,1])
        assert l == [(-1, 2), (0, 3)]

        l = G.get_limits(axes=[1,0])
        assert l == [(0, 3), (-1, 2)]


    def test_meshwidth(self):
        G = self.build_grid()

        m = G.get_meshwidth()
        assert m == array([0.375, 0.3])

        m = G.get_meshwidth(axes=0)
        assert m == array([0.375])

        m = G.get_meshwidth(axes=(1,))
        assert m == array([0.3])

        m = G.get_meshwidth(axes=(0,1))
        assert m == array([0.375, 0.3])

        m = G.get_meshwidth(axes=[1,0])
        assert m == array([0.3, 0.375])


    def test_number_nodes(self):
        G = self.build_grid()

        nn = G.get_number_nodes()
        assert nn == 80

        nn = G.get_number_nodes(axes=0)
        assert nn == 8

        nn = G.get_number_nodes(axes=[1])
        assert nn == 10

        nn = G.get_number_nodes(axes=[0,1])
        assert nn == 80

        nn = G.get_number_nodes(axes=[1,0])
        assert nn == 80

        nn = G.get_number_nodes(overall=False)
        assert nn == array([8, 10])

        nn = G.get_number_nodes(axes=0, overall=False)
        assert nn == array([8])

        nn = G.get_number_nodes(axes=[1], overall=False)
        assert nn == array([10])

        nn = G.get_number_nodes(axes=[0,1], overall=False)
        assert nn == array([8, 10])

        nn = G.get_number_nodes(axes=[1,0], overall=False)
        assert nn == array([10, 8])


    def test_axes(self):
        G = self.build_grid()

        ax = G.get_axes()
        assert len(ax) == 2
        assert ax[0].shape == (8,1)
        assert ax[1].shape == (1,10)

        ax = G.get_axes(axes=0)
        assert len(ax) == 1
        assert ax[0].shape == (8,1)

        ax = G.get_axes(axes=[1])
        assert len(ax) == 1
        assert ax[0].shape == (1,10)


        ax = G.get_axes(axes=[0,1])
        assert len(ax) == 2
        assert ax[0].shape == (8,1)
        assert ax[1].shape == (1,10)

        ax = G.get_axes(axes=[1,0])
        assert len(ax) == 2
        assert ax[0].shape == (1,10)
        assert ax[1].shape == (8,1)

        li = G.get_limits(axes=1)[0]
        ax = G.get_axes(axes=1)[0]
        assert ax.min() == li[0]
        assert ax.max() < li[1]


    def test_grid(self):
        G = self.build_grid()

        no = G.get_nodes()
        assert no.shape == (2, 80)

        li = G.get_limits()
        no = G.get_nodes()
        assert no.min(axis=1)[0] == li[0][0]
        assert no.min(axis=1)[1] == li[1][0]
        assert no.max(axis=1)[0] < li[0][1]
        assert no.max(axis=1)[1] < li[1][1]
