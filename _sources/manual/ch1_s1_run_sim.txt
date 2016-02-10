Set up and run a simulation
---------------------------

Let's first show how to set up a single simulation. The basic work-flow consists
of several steps. First we have to prepare the simulation, then we run the main
simulation program. This gives us a data file with the simulation results. Then
we can apply various post processing steps, for example the computation of
energies, plotting of norms and many more.

The first step is to create a `configuration file` and set the simulation parameters.
We use the file ``examples/harmonic_oscillators/harmonic_1D_p.py`` from the examples
collection. The full content of this file is printed below:

.. include:: ../../../examples/harmonic_oscillators/harmonic_1D_p.py
   :code: python

Now we have to run the main simulation program. This is done by the following command:

::

    Main.py harmonic_1D_p.py

where we have to provide the configuration file as command line option of the ``Main.py``
program. The ``Main.py`` command (like all other commands) supports an online help listing
the available options and switches:

::

    Main.py --help
    usage: Main.py [-h] [-o OUTPUTFILE] [-r [RESULTSPATH]] parametersfile

    positional arguments:
      parametersfile        The simulation configuration parameters file.

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            The data file to write the transformed data.
      -r [RESULTSPATH], --resultspath [RESULTSPATH]
                            Path where to put the results.

When the program terminates, it leaves a file called ``simulation_results.hdf5`` which
contains all the simulation data. This ``hdf`` results file is stored in the `local`
directory where the script was called (unless specified otherwise by the ``-o`` or ``-r``
switches) and `not` in the directory where the configuration file was loaded from.
The script refuses to run if it would overwrite an existing ``hdf`` file.
This measure is in place to prevent you from data loss.

We can use the program ``hdfview`` for example to gain some insight of the contents
of the file. More tools for manipulating `hdf` files can be found on the web-page
of the [hdfgroup]_. Two very useful command line tools are ``h5ls`` and ``h5dump``.

.. [hdfgroup] http://www.hdfgroup.org/products/hdf5_tools/
