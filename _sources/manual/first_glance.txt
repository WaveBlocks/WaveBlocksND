A first glance
==============

Introduction
------------

The *WaveBlocks* project is a collection of reusable software components
providing many of the objects used in the study of semi-classical wavepackets.
Currently it's all about the time-dependent Schroedinger equation and the time
evolution of initial states.

One of the main goals is to provide a set of building blocks - hence the project's
name - that are well tested and reliable. The included features range from
very simple mathematical operations like specialized quadrature rules to basic
data structures for semi-classical wavepackets to more high-level simulation
algorithms and some non-standard plotting functions. Of course there are also
routines included for saving, managing and evaluating simulation results inflexible
manner. All of these components are put together in an easy to use and easy to
extend framework.

The whole project is written in the `Python` programming language and a strong
emphasis is put on readable code and clean software design. Speed and efficiency
are, however, not the main concern for the moment.


Download
--------

The WaveBlocks project has its home at GitHub:

  https://github.com/WaveBlocks

where there are several closely related ``git`` repositories. The latest version
of this project can be found in the ``git`` repository:

  https://github.com/WaveBlocks/WaveBlocksND

Check out the ``master`` branch for the main development version or one of the
tags for stable and released versions. There are also branches used for developing
special features before they are ready to get merged into the ``master`` branch.


Dependencies
------------

The `WaveBlocksND` software depends on some more or less well known `Python` packages,
that are not installed by default. On any recent Linux distribution (for example `Debian`)
you can use the package manager to download and install all the dependencies.
The simulation code is compatible with `Python 3.4` and above. Hence, make sure to
get the correct version of all dependencies. All packages can also be installed
via the ``pip`` `tool <https://pip.pypa.io/en/latest/>`_ from
the `Python Package Index <https://pypi.python.org/pypi>`_.

The necessary dependencies are listed here together with a brief statement of why we need the package.

* ``numpy``, available from https://www.numpy.org

  Numpy provides fast multidimensional arrays.

* ``scipy``, available from https://www.scipy.org

  Scipy interfaces to fast numerical subroutines (BLAS, LAPACK, FFTW).

* ``sympy``, available from https://sympy.org

  Sympy gives rise to symbolic calculation.

* ``matplotlib``, available from http://matplotlib.org

  Matplotlib is used for two-dimensional plotting.

* ``h5py``, available from https://h5py.org

  H5py is the interface to the Hierarchical Data Format `hdf`.

* ``mayavi2``, available from https://github.com/enthought/mayavi

  Mayavi2 is used for three-dimensional plotting.

The ``mayavi2`` dependency is weak in the sense that the whole `WaveBlocks` code
runs fine without it and only the 3-dimensional plotting functions are not available.


Installation
------------

In the following section be sure not to mix up the name `WaveBlocks` that
stands for the whole project and the name ``WaveBlocksND`` which is the
name of the `Python` package.

By using the `Python` ``virtualenv`` `utility <https://virtualenv.pypa.io/en/latest/>`_,
we can easily create a completely isolated `Python` environments. Any `Python`
package can then be installed locally within such an environment without affecting
the outer world.


For Users
~~~~~~~~~

First, we create a new `virtualenv` called ``waveblocks`` which we allow access to
all the `Python` packages that are already installed on the system:

::

    virtualenv -p /usr/bin/python3.4 --system-site-packages waveblocks

This assumes that the dependencies are installed globally on the system.
If we want a truly isolated environment, we can omit the switch ``--system-site-packages``.
This will cause all packages to be built and installed locally to this `virtualenv`.

Next, we switch into the directory and active the `virtualenv`:

::

   cd waveblocks
   source ./bin/activate

Finally, we can install the `WaveBlocks` software directly from the `git` repository
by running:

::

   pip install git+https://github.com/WaveBlocks/WaveBlocksND.git#egg=WaveBlocksND

That is all we have to do. Before using the software, **always** go back to the ``waveblocks`` directory
and run ``source ./bin/activate`` there. Note that this action is **local** to the current shell instance.


For Developers
~~~~~~~~~~~~~~

If you plan to modify the source code, first clone the `git` repository:

::

   git clone https://github.com/WaveBlocks/WaveBlocksND.git ~/waveblockssource

The clone is put into the directory ``~/waveblockssource``. The directory listing should now look like:

::

   LICENSE  MANIFEST.in  README.md  WaveBlocksND  doc  examples  scripts  setup.py

Next, you can create and activate the `virtualenbv` like above. The run the following command
to install the software within the `virtualenv`:

::

   pip install -e ~/waveblockssource

Note the ``-e`` switch, which connects the source code in ``~/waveblockssource`` with the
installation inside the `virtualenv` at ``~/waveblocks``.


Software Overview
-----------------

The `WaveBlocks` package consists of two parts. On one hand there is a `Python` library
and on the other hand there are a bunch of scripts that use this library to implement common
computations.


The Library
~~~~~~~~~~~

The library is called ``WaveBlocksND`` and can be imported in `Python` by the usual procedure:

::

   import WaveBlocksND

this will bring a large number of objects into scope:

::

   dir(WaveBlocksND)
   ['AbstractGrid', 'BasisShape', 'BasisTransformation', 'BasisTransformationHAWP', 'BasisTransformationWF', ...

For the details, please refer to the class documentation.


The Scripts
~~~~~~~~~~~

The scripts (everything in the ``scripts/`` directory) perform simulations, data evaluation and plotting and
can be called from anywhere. These files are just plain `Python` scripts that import the ``WaveBlocksND`` `Python`
module. They all support online help when called with the ``--help`` switch.


Supported platforms
-------------------

The ``WaveBlocksND`` code might run on `Windows` and `OS X` and the various
`BSD` variants too, provided that the required `Python` dependencies are installed.
However, this has never been tested. The primary development platform is `GNU/Linux`
and in particular the `Debian` distribution.
