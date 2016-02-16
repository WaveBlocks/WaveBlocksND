"""The WaveBlocks Project

Provides several computation routines for
handling time and timesteps.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2015, 2016 R. Bourquin
@license: Modified BSD License
"""

from scipy import floor

__all__ = ["TimeManager"]


class TimeManager(object):
    r"""This class performs several computations with time, timesteps and so forth.

    The important quantities here are:

    ============ ============== ======================================================
    Quantity     Parameter Name Description
    ============ ============== ======================================================
    :math:`T`    T              the fixed simulation end time
    :math:`\tau` dt             the size of the timestep
    :math:`N`    nsteps         the overall number of timesteps.
    :math:`t`                   an unspecified time in the interval :math:`[0, T]`
    :math:`n`                   an unspecified timestep in the interval :math:`[0, N]`
    ============ ============== ======================================================

    The important relations that hold are :math:`T = N \tau` and
    in analogy :math:`t = n \tau`. There are also conversion routines
    for :math:`t` and :math:`n`.

    The simulation parameters handed over to the constructor must contain at least
    two out of the three values :math:`T`, :math:`\tau` and :math:`N`. If all three
    are given, the user is responsible for compatible values.

    Additionally the class contains some routines for determining
    if and when some events (for example saving data) should occur.
    """

    def __init__(self, parameters):
        if parameters is None:
            parameters = {}

        # We need two out of three: T, dt and nsteps
        have_enough = 0

        if "T" in parameters:
            self._T = float(parameters["T"])
            have_enough += 1
        else:
            self._T = None

        if "dt" in parameters:
            self._dt = float(parameters["dt"])
            have_enough += 1
        else:
            self._dt = None

        if "nsteps" in parameters:
            self._nsteps = int(parameters["nsteps"])
            have_enough += 1
        else:
            self._nsteps = None

        if have_enough < 2:
            raise KeyError("Parameters provide to little data to construct a 'TimeManager'.")

        if self._T is None:
            self._T = self.compute_endtime()

        if self._dt is None:
            self._dt = self.compute_timestep_size()

        if self._nsteps is None:
            self._nsteps = self.compute_number_timesteps()

        # Interval for regular events
        self._interval = 1
        if "write_nth" in parameters:
            self.set_interval(int(parameters["write_nth"]))

        # List of timesteps of irregular events
        self._eventtimes = []
        if "save_at" in parameters:
            self.add_to_eventlist(parameters["save_at"])


    def __str__(self):
        s  = "TimeManager configured with:\n"
        s += " Final time     T: "+str(self._T)+"\n"
        s += " Timestep size dt: "+str(self._dt)+"\n"
        s += " Number of steps : "+str(self._nsteps)+"\n"
        return s


    def set_T(self, T):
        r"""Set the simulation endtime :math:`T`.

        :param T: The simulation end time.
        """
        self._T = float(T)


    def set_dt(self, dt):
        r"""Set the simulation timestep size :math:`\tau`.

        :param dt: The simulation timestep size.
        """
        self._dt = float(dt)


    def set_nsteps(self, nsteps):
        r"""Set the number of timesteps the simulation runs.

        :param nsteps: The number :math:`n` timesteps we do.
        """
        self._nsteps = int(nsteps)


    def get_T(self):
        r"""Set the simulation endtime :math:`T`.

        :returns: The endtime :math:`T`.
        """
        return self._T


    def get_dt(self):
        r"""Get the simulation timestep size :math:`\tau`.

        :returns: The timestep :math:`\tau`.
        """
        return self._dt


    def get_nsteps(self):
        r"""Get the number :math:`n` of timesteps the simulation runs.

        :returns: the number :math:`n` of timesteps.
        """
        return self._nsteps


    def compute_endtime(self):
        r"""Computes the simulation endtime :math:`T`.

        :returns: The endtime :math:`T`.
        """
        if self._T is not None:
            return self._T
        else:
            return float(self._nsteps * self._dt)


    def compute_timestep_size(self):
        r"""Computes the simulation timestep size :math:`\tau`.

        :returns: The timestep :math:`\tau`.
        """
        if self._dt is not None:
            return self._dt
        else:
            return self._T / float(self._nsteps)


    def compute_number_timesteps(self):
        r"""Computes the number :math:`n` of time steps we will perform.

        :returns: the number :math:`n` of timesteps.
        """
        if self._nsteps is not None:
            return self._nsteps
        else:
            return int(floor(self._T / float(self._dt)))


    def compute_timestep(self, t):
        r"""Compute the timestep :math:`n` from a time :math:`t` such that
        :math:`t = n \tau` holds.

        :param t: The time t of which we want to find the timestep number.
        :returns: The corresponding timestep :math:`n`.

        Note that the user has to ensure that time :math:`t` is an integral
        multiple of :math:`\tau`.
        """
        stepo = t / self._dt
        step = round(stepo)

        if abs(stepo - step) > 1e-10:
            print("Warning: Questionable rounding for timestep computation!")

        return int(step)


    def compute_time(self, n):
        r"""Compute the time :math:`t` from a timestep :math:`n` such that
        :math:`t = n \tau` holds.

        :param n: The timestep n of which we want to find the corresponding time.
        :returns: The corresponding time :math:`t`.
        """
        return float(n * self._dt)


    def set_interval(self, interval):
        r"""Set the interval for regular events.

        :param interval: The interval at which regular events get triggered.

        Note that a value of ``0`` means there are no regular events.
        """
        self._interval = int(interval)


    def add_to_eventlist(self, alist):
        r"""Add a list of times and/or timesteps to the list of
        times when irregular events get triggered.

        :param alist: A list with integers (interpreted as timesteps)
                      and/or floats (interpreted as times)

        Note that the times and timesteps can be mixed and need not to be
        given in monotone order.
        """
        timesteps = []

        # If the list is empty (global default), shortcut
        if len(alist) == 0:
            return

        # Integers are interpreted as timesteps, floats are interpreted as times (and converted to timesteps)
        for item in alist:
            if type(item) == int:
                timesteps.append(item)
            elif type(item) == float:
                timesteps.append(self.compute_timestep(item))

        # Validate timesteps and check if n in [0,...,N]
        tmp = len(timesteps)
        timesteps = [i for i in timesteps if 0 <= i <= self._nsteps]

        if tmp != len(timesteps):
            print("Warning: Dropped %d timestep(s) due to invalidity!" % (tmp - len(timesteps)))

        # Assure unique elements, just silently remove duplicates
        oldlist = set(self._eventtimes)
        newlist = set(timesteps)
        times = list(oldlist.union(newlist))
        # Sort in ascending order
        times.sort()
        # Write back
        self._eventtimes = times


    def compute_number_events(self):
        r"""Compute the number of events we will perform during the simulation.
        This can for example be used to determine how much space to allocate
        in the output files if the events are times at which simulation data
        is saved.

        :returns: The number of events.
        """
        # We do not save at regular intervals
        if self._interval == 0:
            # The number of saves resulting from saving at a regular interval is zero
            n_si = 0
            # Determine the number of saves resulting from the savelist
            n_sl = len(self._eventtimes)
        # We do save at regular intervals
        else:
            # Determine the number of saves resulting from saving at a regular interval
            n_si = 1 + self._nsteps // self._interval
            # Determine the number of saves resulting from the savelist and
            # exclude the timesteps which coincide with the regular intervals
            n_sl = len([i for i in self._eventtimes if i % self._interval != 0])

        # Total number of saves we will perform is given by the sum
        number_events = n_si + n_sl

        return number_events


    def is_event(self, n):
        r"""Determine if an event occurs right now.

        :param n: The current timestep in question.
        :returns: ``True`` or ``False``.
        """
        if self._interval == 1:
            # Save every timestep
            return True
        elif self._interval != 0 and n % self._interval == 0:
            # Save every k-th timestep specified by the interval
            return True
        elif n in self._eventtimes:
            # Save if the n is in the list of timesteps
            return True

        return False
