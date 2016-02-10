Running multiple simulations
----------------------------

Now we know how to run a single simulation. But most of the time we want
to run a multitude of simulations. This is not more difficult, only the work-flow
changes a little bit. Throughout the next section we work in an arbitrary
directory. All files referenced are assumed to lie within this working directory.

Preparation and Meta-configurations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we need to generate a bunch of configurations. Of course we could write
all the files by hand. However, for a set of simulations where just one or a
few parameters vary, we can avoid this tedious work. The tool that takes over
the task is named ``ConfigurationGenerator.py``. It takes a so called `meta configuration`
and then produces a set of ordinary configuration files. The synopsis for this
tool is:

::

    ConfigurationGenerator.py --help
    usage: ConfigurationGenerator.py [-h] [-d DESTINATION] metaconfiguration

    positional arguments:
      metaconfiguration     The meta-configuration file.

    optional arguments:
      -h, --help            show this help message and exit
      -d DESTINATION, --destination DESTINATION
                            The destination where to store the configurations
                            generated.

Let's look at a simple example: assume that our sample meta configuration file
is ``metaconfiguration_02.py``, its content is reprinted below:

::

    # Global parameters that stay the same for all simulations :
    GP = {}
    GP["algorithm"] = "\"fourier\""
    GP["potential"] = "\"delta_gap\""
    GP["T"] = 3
    GP["dt"] = 0.02
    GP["parameters"] = "[ (1.0j, 1.0-6.0j, 0.0, 1.0, -6.0), (1.0j, 1.0-6.0j, 0.0, 1.0, -6.0) ]"
    GP["coefficients"] = [ [(0 ,1.0)], [(0,0.0)] ]
    GP["basis_size"] = 2
    GP["ngn"] = 2**12
    GP["f"] = 4.0
    GP["write_nth"] = 2

    # Local parameters that change with each simulation
    LP = {}
    LP["eps"] = [0.1, 0.5]
    LP["delta"] = ["0.5*eps", "1.0*eps", "1.5*eps"]

The file is just another plain `Python` file with only informal constraints.
There must be two dicts named ``GP`` and ``LP`` in the top level name-space.
The first one, ``GP``, contains all the parameters that are `global` to the
set of configuration. While the second one, ``LP``, contains lists of the
parameters that vary with each simulation. The configuration generator then
computes the Cartesian product of all these lists in ``LP``. Then, for each
tuple of this Cartesian product it adds all parameters from ``GP``, yielding
a single configuration. Additionally to these two variables there can be
another one which is used for global preambles. This variable has to be called
``PA`` and holds a (multi-line) `Python` string of valid `Python` code. These
statements are written to the very top of every configuration file generated.

We can run the configuration generator as:

::

    ConfigurationGenerator.py metaconfiguration_02.py

and it will create the directory ``autogen_configurations`` where it puts
all the configuration files. Let's take a look into this directory:

::

    ls -l autogen_configurations/

prints:

::

    Parameters[eps=0.1][delta=0.5eps].py
    Parameters[eps=0.1][delta=1.0eps].py
    Parameters[eps=0.1][delta=1.5eps].py
    Parameters[eps=0.5][delta=0.5eps].py
    Parameters[eps=0.5][delta=1.0eps].py
    Parameters[eps=0.5][delta=1.5eps].py

and we find 6 configuration files. One file for each combination of a value for
``eps`` and one for ``delta``. The file names contain all local parameters as ``key=value``
pairs. These can be used later in the post processing step by the functions from
the :py:class:`FileTools` sub-module for sorting and grouping the simulations with
respect to almost arbitrary criteria.

These configuration files can now be fed to the main simulation program one
after another as shown in the last section. We could again do this manually but
there is a better solution.


The batch loop
~~~~~~~~~~~~~~

There is a simple `Python` script called ``BatchLoop.py`` which does nothing else than running
simulations for a set of configurations. The usage is really simple.

::

    BatchLoop.py --help
    usage: BatchLoop.py [-h] -c CONFIGURATIONS [-r RESULTSPATH] [-m MAXWORKERS]

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIGURATIONS, --configurations CONFIGURATIONS
                            Path to the 'configuration' directory.
      -r RESULTSPATH, --resultspath RESULTSPATH
                            Path to the 'results' directory.
      -m MAXWORKERS, --maxworkers MAXWORKERS
                            Maximal number of parallel jobs.

We can run as many simulations as we like. Each simulation is run independently
from all others and there is a limit of ``MAXWORKERS`` simulations run in parallel.

.. warning:: The ``BatchLoop.py`` command runs only with `Python` 3.4 or later because
             of the use of new library features providing support for concurrent execution
             of code!

We have to provide a directory where the results should end up:

::

    mkdir results

Now it is time to call the ``BatchLoop.py`` script. The simple call looks like:

::

    BatchLoop.py -c autogen_configurations -r results

This will create new directories in ``results`` whose names correspond to the
configuration files used. It will call the ``Main.py`` script for each simulation
configuration provided. After this it will run a bunch of data computation and plotting
scripts. If we now look into the results directory by:

::

    ls results

we see the listing:

::

    Parameters[eps=0.1][delta=0.5eps]
    Parameters[eps=0.1][delta=1.0eps]
    Parameters[eps=0.1][delta=1.5eps]
    Parameters[eps=0.5][delta=0.5eps]
    Parameters[eps=0.5][delta=1.0eps]
    Parameters[eps=0.5][delta=1.5eps]

and for the results of a single simulation (notice the necessary shell character
escapes, you can also write the name without escapes in a pair of ``"``.)

::

    ls results/Parameters\[eps\=0.1\]\[delta\=0.5eps\]

we have the following bunch of files:

::

    energies_block0.png
    energy_drift_block0.png
    norms_block0.png
    norms_drift_block0.png
    norms_sqr_block0.png
    Parameters[eps=0.1][delta=0.5eps].py
    simulation_results.hdf5

Each directory within results contains at least the simulation parameters
file (``Parameters[eps=0.1][delta=0.5eps].py``) and the simulation results
file (``simulation results.hdf5``). If there were some plots generated,
then these files are here too.


Running more scripts
~~~~~~~~~~~~~~~~~~~~

Sometimes you may wish to run a script for a set of simulations long after the
batch loop has terminated. Maybe you decided to compute a new observable or
whatever. It would be tedious to call the script with each ``simulation_results.hdf5``
and its correct file path manually. Exactly for this reason there is a script named
``ForAll.py``. For example assume we want to plot the potential used in each simulation
(which is identical in our example but never mind). Then we call:

::

    ForAll.py PlotPotential.py

which starts by printing:

::

    Will execute the code in 'PlotPotential.py' for all files in 'results'
     Executing code for datafile in results/Parameters[eps=0.5][delta=1.0eps]
     ...

and after a while quits with the text ``Done`` on the last output line. The script
can take the path of the directory where the results lie (in the example above
this is ``./results/``) as a third command line argument.
