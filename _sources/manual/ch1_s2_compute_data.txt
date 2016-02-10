Computing more data
-------------------

After we have run a simulation the output file ``simulation_results.hdf5``
contains all data that were computed during the simulation. This is for example
wave-function values or wave-packet parameters etc. depending on the exact setup
run. Usually we want also to compute some properties of the time evolution. This
is done in a second step called `post processing` of the data. There are several
scripts in the ``scripts/`` sub-directory which post-process the simulation data.

Assume we want to compute the norms and energies of the wave function during its
time evolution. These properties are not computed while running the simulation,
but we can get them easily from the stored information. The following sections
will show how to compute these data and store them in the output file
``simulation_results.hdf5`` too.

All post-processing and plotting scripts can be called with an argument ``--help``
and provide modern command line switch handling.

::

    ComputeNorms.py --help

and will print a help message:

::

    usage: ComputeNorms.py [-h] [-d [DATAFILE]] [-b [BLOCKID [BLOCKID ...]]]
                           [-r [RESULTSPATH]] [-et]

    optional arguments:
      -h, --help            show this help message and exit
      -d [DATAFILE], --datafile [DATAFILE]
                            The simulation data file.
      -b [BLOCKID [BLOCKID ...]], --blockid [BLOCKID [BLOCKID ...]]
                            The data block to handle.
      -r [RESULTSPATH], --resultspath [RESULTSPATH]
                            Path where to put the results.
      -et, --eigentransform
                            Transform the data into the eigenbasis before
                            computing norms.

Norms
~~~~~

Computing norms is trivial and fast. Just run the script:

::

    ComputeNorms.py

This will compute the norms of all wave-packets or wave functions
depending on what the simulation setup was and what is already stored
in ``simulation_results.hdf5``.

Energies
~~~~~~~~

Asking for the energies is almost equally trivial as computing norms.
All we need is to run:

::

    ComputeEnergies.py

which will compute kinetic and potential energies.

Autocorrelations
~~~~~~~~~~~~~~~~

The computation of auto-correlations is a bit more complicated. What
we want to compute is the following overlap integral (here discussed
in case of wave-packets):

.. math::
   \langle \Psi(0) | \Psi(t) \rangle

which compares the wave-packet at time :math:`t` with the initial value
:math:`\Psi(0)` at time 0. Because this involves wave-packets at two different
times we need a specialized quadrature to get accurate results.  We have to tell
the script which quadrature we would like to use. This is done best by adding a
top-level snippet like the following to the original simulation setup
configuration *before* the simulation is run. This will choose the
:py:class:`NSDInhomogeneous` quadrature transformation using
:py:class:`GaussHermiteOriginalQR` with 4 nodes and is for a one-dimensional
setup:

::

    # Configurations needed for computation of observables
    observables = {
        "autocorrelation" : {
            "innerproduct" : {
                "type" : "InhomogeneousInnerProduct",
                "delegate" : {
                    "type" : "NSDInhomogeneous",
		    "qr" : {"dimension": 1, "order": 4, "type": "GaussHermiteOriginalQR"}
                    }
                }
            }
        }

.. warning::
   It is essential to take :py:class:`GaussHermiteOriginalQR` as quadrature rule
   used by the :py:class:`NSDInhomogeneous` transformation.

As a second example we show the corresponding snippet in case of a three
dimensional simulation setup:

::

    # Configurations needed for computation of observables
    observables = {
        "autocorrelation" : {
            "innerproduct" : {
                "type" : "InhomogeneousInnerProduct",
                "delegate" : {
                    "type" : "NSDInhomogeneous",
                    "qr": {
                        "type": "TensorProductQR",
                        "dimension": 3,
                        "qr_rules": [
                            {"dimension": 1, "order": 5, "type": "GaussHermiteOriginalQR"},
                            {"dimension": 1, "order": 5, "type": "GaussHermiteOriginalQR"},
                            {"dimension": 1, "order": 5, "type": "GaussHermiteOriginalQR"}],
                        }
                    }
                }
            }
        }

The only thing we have to do then is to call the corresponding post-processor script:

::

    ComputeAutocorrelation.py


Wave-packet sampling
~~~~~~~~~~~~~~~~~~~~

If we made a simulation with wave-packets only and want to sample them
on a regular grid for example for plotting then there is a script for this purpose:

::

    usage: ComputeEvaluateWavepacketsCanonical.py [-h] [-d [DATAFILE]]
                                                  [-b [BLOCKID [BLOCKID ...]]]
                                                  [-p [PARAMETERSFILE]]
                                                  [-r [RESULTSPATH]] [-et]

    optional arguments:
      -h, --help            show this help message and exit
      -d [DATAFILE], --datafile [DATAFILE]
                            The simulation data file.
      -b [BLOCKID [BLOCKID ...]], --blockid [BLOCKID [BLOCKID ...]]
                            The data block to handle.
      -p [PARAMETERSFILE], --parametersfile [PARAMETERSFILE]
                            The configuration parameter file.
      -r [RESULTSPATH], --resultspath [RESULTSPATH]
                            Path where to put the results.
      -et, --eigentransform
                            Transform the data into the eigenbasis before
                            computing norms.

Eigentransformations
~~~~~~~~~~~~~~~~~~~~

For potentials with multiple energy levels it matters in which basis we compute
observables. Since the simulation is done in the canonical basis and the
observables usually should be computed in the eigenbasis there is a
transformation involved. The scripts shown above do this transformation
internally and there is no need to worry.

However, in case we explicitly do not want the transformation to take place
(for example when working with single-level potentials) there are suitable
post-processing scripts which can be recognized by a ``NET`` in their name:

::

    ComputeNormsNET.py
    ComputeEnergiesNET.py
    ComputeAutocorrelationNET.py

The ``NET`` (No-Eigen-Transformation) variants never do a basis transformation
and compute the requested observables on the data given assuming a correct
basis. There is also a ``CAN`` variant which computes explicitly in the
canonical basis:

::

    ComputeEnergiesCAN.py

The reason why this script exists is that it makes a difference whether
we use :math:`V(x)` or :math:`\Lambda(x)` in the code.


Explicit Eigentransformation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In case we want to convert all the simulation data (think: wave-function values
or wave-packet data) once to the eigenbasis there is this script:

::

    ComputeTransformToEigen.py --help

According to its help text:

::

    usage: ComputeTransformToEigen.py [-h] [-i INPUTFILE] [-o OUTPUTFILE]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUTFILE, --inputfile INPUTFILE
                            The data file to read the data from.
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            The data file to write the transformed data.

it will read the input file ``simulation_results.hdf5`` and write output into a
new data file. A typical invoke could look like:

::

    ComputeTransformToEigen.py -i simulation_results.hdf5 -o simulation_results_eigen.hdf5
