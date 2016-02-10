Data storage
============

What data are written to disk. How can we retrieve data, IOM basics, usage, etc


How IOM works
-------------

The so-called `IOManager` is responsible for storing all our data. It provides a
meaningful API for storing and retrieving simulation data and the goal is to
make data handling from scripts as easy as possible. The IOManager uses the low-level
``hdf5`` file format to actually store the numerical data efficiently. Dealing directly
with the hdf5 API provided by ``h5py`` would be cumbersome as we would have
to remember much more details about how the data are stored inside an hdf file.
With this thin layer we just tell the IOM which data we want to store or load and
it performs all the low-level stuff behind our back.

Please note that the tab-completion of ``ipython`` won't work as usual
on ``IOManager`` instances because of its plugin architecture. The plugins
allow to add functionality at runtime and only when its really used. Thus a
(member)function may be loaded right at the moment it gets called the first time.
This is the reason why tab-completion and introspection will not work for
(member)functions that had never been called before.


What gets stored
----------------

Each file containing simulation results is basically divided into `datablocks`.
There is one special block called the `global datablock` which stores
data that are identical for the whole simulation (for example space domain grids,
simulation parameters etc). Then there can be an arbitrary number of normal data
blocks which can store various data related to wavepackets, wavefunctions and observables.
Each of these data sets is optional and there are functions to query if specified
data is available. The next figure shows the coarse structure of any simulation
results file.

.. figure:: ./fig/structure_result_file.png
   :align: center

   Coarse structure of a file containing simulation results.

The figure below shows the internal structure of a single data block.
In this structure not all data objects always exist depending on what
computations were performed. The dark blocks are at the level of individual data
tensors while the lighter grey boxes represent hdf groups. Note that not all data
sets may exist at all and that each group can have different subsets. For example
if you never computed observables, then this entire block is missing. The
wavefunction data can come from a simulation with the Fourier propagator or from
the evaluation of wavepackets on a given domain-wide grid.

.. figure:: ./fig/structure_datablock.png
   :align: center

   Possible structure of a single data block. Not all data always exist.


Saving data at times and timesteps
----------------------------------

Storing simulation data can happen in various different ways. For example you
can store data at regular time intervals. Or at a list of fixed points in time.
Both is easily possible with the tools provided by the ``IOManager`` together
with the ``TimeManager``. While the ``IOManager`` is responsible for
saving and loading the data, the ``TimeManager`` is used for all computations
related with time, timesteps and so on, for example to convert a list of times
into a list of timesteps or checking if a given time is is within the simulated
time range etc.

The two parameters ``write_nth`` and ``save_at`` are used to configure the
way you wish to save data. While the first is used to specify the details of saving
at regular time intervals, the second one provides the means to specify a list
of points in time when saving should take place. A few examples of saving at regular
intervals::

  # Save data at each timestep
  write_nth = 1

  # Save data each 5th timestep
  write_nth = 5

  # Never save data
  write_nth = 0

Please note that this scheme is rigid in the sense that if for example the timestep
corresponding to the end of the simulation is not an integer multiple of the value
of this parameter then the data from the end is missing. (This should be quite obvious!)

The parameter ``save_at`` has to be a python list containing integers
and/or floats. There is a *big difference* between the two data types
you always have to be aware of! Integer values are interpreted as `timesteps`
while floats will be taken as `times`. A few examples on saving at specified
times only::

  # Save at timestep 3, 6, 7, 13 and 19
  save_at = [3, 6, 7, 13, 19]

  # Save at the end time only
  # Assuming T = 5.34 and T is an integer multiple of dt!
  save_at = [5.34]

  # Save at a few times
  # This is usefull to compare simulation results of simulations
  # with different timestep sizes. Of course the times have to be
  # integer multiples of *all* timestep sizes in consideration!
  save_at = [3.2, 4.5, 8.7, 19.3]

You can freely mix the two approaches and specify crazy things like
the following::

  write_nth = 15
  save_at = [1, 2, 3, 4.5, 10, 3.2, 40, 23.45, 23.55]

which translates to the: `Save the data each 15 steps and additionally
save the data at the timesteps 1, 2, 3, 10 and 40 and save the data at the time 3.2,
23.45 and 23.55.` It is assumed that `time` is an integer multiple of the
``timestep`` size. (Otherwise more or less careful rounding will be applied.)
The list doesn't have to be in monotone order and duplicates will be removed as well
as values outside the interval :math:`[0, T]` where :math:`T` is the time at which
the simulation stops. A good use case for a mixed specification is for example saving at big
intervals but including the very end of the simulation::

  write_nth = 35
  save_at = [5.34]    # Same assumption as above

Note that even if you disable saving data entirely be setting::

  write_nth = 0     # Default is 1
  save_at = []      # Default is []

you will end up with a hdf5 file still containing the initial values as they
are at time equal 0 (before the first timestep was made).


Retrieving the simulation parameters
------------------------------------

From a hdf5 file with the simulation data we can get back the parameters this
simulation used. Retrieval is trivial, the following commented interactive python
session shows the basics which can of course be used in a user script too::

  >>> from WaveBlocks import IOManager
  >>> iom = IOManager()                         # create an IOM instance
  >>> iom.load_file("simulation_results.hdf5")  # load the data file
  >>> sim_params = iom.get_parameters()         # request the parameters
  >>> print(sim_params)
  ====================================
  Parameters of the current simulation
  ------------------------------------
  [...]

With only three trivial lines of code we get back all the parameters
that were used for the simulation!


Load simulation data
--------------------

Simulation data can be loaded from a given ``simulation_results.hdf5`` file by
an IOManager instance. You can even do this inside an interactive ``ipython``
session. The API is quite trivial, all functions for loading data have their name
prefixed by ``load_`` as for example in ``load_energy(...)``. Every function
for loading and saving data has a keyword argument ``block`` defaulting to 0
which tells the IOManager from which data block to take the requested data.
For quantities that represent time series, the load functions also provide a keyword
argument ``timestep`` that can be used to load data from a single timestep.
The default is ``None`` meaning `load the data from all timesteps`.
A sample of such an interactive session could look like this::

  >>> from WaveBlocks import IOManager
  >>> iom = IOManager()                          # Create a new IOManager instance
  >>> iom.open_file("simulation_results.hdf5")   # And open a given hdf5 file

  >>> print(iom)
    IOManager instance with open file simulation_results.hdf5

  >>> ekin, epot = iom.load_energy()         # Load the energies from a simulation
    Requested function: load_energy          # Don't bother about the messages
    Plugin to load: IOM_plugin_energy        # concerning the plugins.

  >>> ekin.shape                             # We see the the energies are given
    (301, 1)                                 # as time series over 301 timesteps
  >>> epot.shape
    (301, 1)

  >>> tg = iom.load_energy_timegrid()        # Load the corresponding timegrid which
                                             # contains the timesteps when the data
  >>> tg.shape                               # was saved. This is important if the
    (301,)                                   # data was saved at non-regular intervalls.

  >>> iom.finalize()                         # Close the hdf5 file

  >>> plot(tg, ekin)                         # Plot the kinetic energy over time

Of course all this works exactly the same inside any regular python script.
For a complete list of all the ``load_`` functions please see the API
documentation or the docstrings.


Working with simulation data
----------------------------

The following code snippet shows how to perform a data transformation task
for all blocks of a simulation results file.::

  >>> iom = IOManager()
  >>> iom.open_file("testdata.hdf5")

  >>> for blockid in iom.get_block_ids():      # Iterate over all data blocks
          if iom.has_energy(block=blockid):    # If the current data block containes
              ...                              # energies we may do something
