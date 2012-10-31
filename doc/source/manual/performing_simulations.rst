Using `WaveBlocks` for performing simulations
=============================================

In this chapter we show how to use the `WaveBlocks` framework for performing
simulations. The process is always the same and consists of a pre-processing
step, a main step and a post-processing step. The preprocessing step is where
we configure the simulations we want to perform. Then there is the main step
where the simulations are run. Finally, there follows a postprocessing step where
we evaluate the data and (optionally) create visualisations. We will see that the
post processing step consists of many small and independent substeps reflecting
the various options of what to do with the data obtained.

Set up and run a single simulation
----------------------------------

Let's first show how to set up a single simulation. The basic workflow consists
of several steps. First we have to prepare the simulation, then we run the main
simulation program. This gives us a data file with the simulation results. Then
we can apply various post processing steps, for example the computation of
energies, plotting of norms and many more.

The first step is to create a `configuration file` and set the parameters. Let's call
the file ``parameters_01.py``. The full content of this file is printed below

::

    algorithm = "hagedorn"

    T = 12
    dt = 0.01

    dimension = 1
    ncomponents = 1

    eps = 0.1

    potential = "quadratic"

    # The parameter set of the initial wavepacket
    Q = [[1.0]]
    P = [[1.0j]]
    q = [[1.0]]
    p = [[0.0]]
    S = [[0.0]]

    # What it takes to specify a wavepacket!
    wp0 = {
        "type" : "HagedornWavepacket",
        "dimension" : 1,
        "ncomponents": 1,
        "eps" : 0.1,
        "Pi" : [q,p,Q,P,S],
        "basis_shapes" : [{
            "type" : "HyperbolicCutShape",
            "K" : 10,
            "dimension" : 1
        }],
        "coefficients" : [[ ((0,), 1.0) ]],
        "quadrature" : {
            "type" : "HomogeneousQuadrature",
    	'qr': {
                'type': 'TensorProductQR',
                'dimension': 1,
                'qr_rules': [{'dimension': 1, 'order': 14, 'type': 'GaussHermiteQR'}]
            }
        }
    }

    # Which wavepackets are initial values
    initvals = [ wp0 ]

    leading_component = 0

    # How often do we write data to disk
    write_nth = 5

    matrix_exponential = "pade"

For an overview of the available settings, see ??.

Now we have to run the main simulation program. This is done by the following
command

::

    python Main.py parameters_01.py

where we have to provide the configuration file as the first command line option
of the ``Main.py`` program. When the program terminates, it leaves a file called
``simulation_results.hdf5`` which contains all the simulation data. We can use
the program ``hdfview`` to gain some insight of the contents of the file.

Now we can start with the post processing of the data. Assume we want to plot
the norms and energies of the wave function during the time evolution. These
are not computed during the simulation, but we can get them from the stored
information. The following two commands will compute these data and store
them in ``simulation_results.hdf5``

::

    python ComputeNorms.py
    python ComputeEnergies.py

What remains is the visualization the data. This is done by two plot scripts

::

    python PlotNorms.py
    python PlotEnergies.py

The post processing step usually splits into two substeps. First we compute
additional data and then we visualise these data. The two substeps are performed
by individual scripts. All these scripts optionally take the filename or filepath
of the ``simulation_results.hdf5`` as a further command line argument.


Running multiple simulations
----------------------------

Now we know how to run a single simulation. But most of the time we want
to run a multitude of simulations. This is not more difficult, only the workflow
changes a little bit. Throughout the next section we work in an arbitrary
directory. All scripts called and all files referenced are assumed to lie within this
working directory.

Preparation and Meta-configurations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we need to generate a bunch of configurations. Of course we could write
all the files by hand. However, for a set of simulations where just one or a
few parameters vary, we can avoid this tedious work. The tool that takes over
the task is named ``ConfigurationGenerator.py``. It takes a so called `meta
configuration` and then produces a set of ordinary configuration files.

Let's look at a simple example: assume that our sample meta configuration file
is ``metaconfiguration_02.py``, its content is reprinted below

::

    # Global parameters that stay the same for all simulations :
    GP = { }
    GP [ " algorithm " ] = " \" fourier \" "
    GP [ " potential " ] = " \" delta_gap \" "
    GP [ " T" ] = 3
    GP [ " dt " ] = 0 . 02
    GP [ " parameters " ] = " [ (), () ] "
    GP [ " coefficients " ] = [ [ (0 , 1 . 0 ) ] , [ (0 , 0 . 0 ) ]
    GP [ " basis_size " ] = 2
    GP [ " ngn " ] = 2 * * 12
    GP [ " f" ] = 4 . 0
    GP [ " write_nth " ] = 2

    # Local parameters that change with each simulation
    LP = { }
    LP["eps"] = [0.1, 0.5]
    LP["delta"] = ["0.5*eps", "1.0*eps", "1.5*eps"]

The file is just another plain python file with only informal constraints.
There must be two dicts named ``GP`` and ``LP`` in the top level namespace.
The first one, ``GP``, contains all the parameters that are `global` to the
set of configuration. While the second one, ``LP``, contains lists of the
parameters that vary with each simulation. The configuration generator then
computes the cartesian product of all these lists in ``LP``. Then, for each
tuple of this cartesian product it adds all parameters from ``GP``, yielding
a single configuration. We can run the configuration generator as:

::

    python ConfigurationGenerator.py metaconfiguration_02.py

and it will create the directory ``autogen_configurations`` where it puts
all the configuration files. Let's take a look into this directory:

::

    ls -l autogen_configurations/

prints

::

    Parameters[eps=0.1][delta=0.5eps].py
    Parameters[eps=0.1][delta=1.0eps].py
    Parameters[eps=0.1][delta=1.5eps].py
    Parameters[eps=0.5][delta=0.5eps].py
    Parameters[eps=0.5][delta=1.0eps].py
    Parameters[eps=0.5][delta=1.5eps].py
