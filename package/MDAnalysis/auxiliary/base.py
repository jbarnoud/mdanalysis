# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Auxiliary Readers --- :mod:`MDAnalysis.auxiliary.base`
======================================================

Base classes for deriving all auxiliary data readers.

.. autoclass:: AuxReader
   :members:

.. autoclass:: AuxFileReader
   :members:

"""

import six

import numpy as np
import math

from ..lib.util import asiterable

from . import _AUXREADERS


class _AuxReaderMeta(type):
    # auto register on class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for f in fmt:
                _AUXREADERS[f] = cls


class AuxReader(six.with_metaclass(_AuxReaderMeta)):
    """ Base class for auxiliary readers.

    Allows iteration over a set of data from a trajectory, additional 
    ('auxiliary') to the regular positions/velocities/etc. This auxiliary 
    data may be stored in e.g. an array or a separate file.

    Auxiliary data may be added to a trajectory by 
    :meth:`MDAnalysis.coordinates.base.Reader.add_auxiliary`, passing either an 
    AuxReader instance or the data/filename, in which case an appropriate reader 
    will be selected based on type/file extension.
   
    The AuxReader will handle alignment of the auxiliary and trajectory 
    timesteps (auxiliary data steps are 'assigned' to the closest trajectory 
    timestep), and for each trajectory timestep will provide a 'representative'
    auxiliary value (or values) based on the steps assigned to that timestep.


    Paramaters
    ----------
    name : str
        Name for auxiliary data. When added to a trajectory, the representative 
        auxiliary value(s) for the timestep are stored as ``ts.aux.name``.
    represent_ts_as : {'closest', 'average'}
        How to calculated representative value of auxiliary data for a 
        trajectory timestep. Currently available:
          *'closest': value from step closest to the trajectory timestep
          *'average': average of values from auxiliary steps assigned to 
                      the trajectory timestep.
    cutoff : float, optional
        Auxiliary steps further from the trajectory timestep than *cutoff* 
        will be ignored when calculating representative values (the default
        value is -1, which indicates all auxiliary steps assigned to that 
        timestep will be used).
    dt : float, optional
        Change in time between auxiliary steps (in ps). If not specified, will
        attempt to be determined from auxiliary data; otherwise defaults to 1ps.
    initial_time : float, optional 
        Time of first auxilairy step (in ps). If not specified, will attempt to
        be determined from auxiliary data; otherwise defaults to 0ps.
    time_col : int, optional
        Index of column in auxiliary data storing time (default value ``None``).
    data_cols : list of str, optional
        Indices of columns containing data of interest to be stored in 
        `step_data` (defaults to all columns).
    constant_dt : bool, optional
        If true, will use dt/initial_time to calculate time even when time
        stored in auxiliary data (default value is ``True``).



    Attributes
    ----------
    step : int
        Number of the current auxiliary step, starting at 0.
    n_steps : int
        Total number of auxiliary steps
    n_cols : int
        Number of columns of data for each auxiliary step.
    time : float
        Time of current auxiliary step, as read from data (if present) or 
        calculated using `dt` and `initial_time`.
    times : list of float
        List of the times of each auxiliary step.
    step_data : ndarray
        Value(s) from the auxiliary data column(s) of interest (per `data_cols`) 
        for the current step.
      
    ts_data : ndarray
        List of 'step_data' from each auxiliary step assigned to the current
        trajectory timestep.
    ts_rep : list of float
        Represenatative value of auxiliary data for current trajectory timestep.

    """
      
    def __init__(self, auxname, represent_ts_as='closest', cutoff=-1, 
                 dt=None, initial_time=None, time_col=None, data_cols=None, 
                 constant_dt=True):

        self.name = auxname

        # update when add new options
        represent_options = ['closest', 'average']
        if represent_ts_as not in represent_options:
            raise ValueError("{0} is not a valid option for calculating "
                             "representative value(s). Enabled options are: "
                             "{1}".format(represent_ts_as, represent_options))
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff

        # set initially to avoid error on first read
        self.n_cols = None

        self.ts_data = None
        self.ts_rep = None

        self._initial_time = initial_time
        self._dt = dt
        self.constant_dt = constant_dt

        self.step = -1
        self._read_next_step()
        self.n_cols = len(self._data)

        if time_col >= self.n_cols:
            raise ValueError("Index {0} for time column out of range (num. "
                             "cols is {1})".format(time_col, self.n_cols))
        if data_cols:
            for col in data_cols:
                if col >= self.n_cols:
                    raise ValueError("Index {0} for data column out of range ("
                                     "num. cols is {1})".format(col,self.n_cols))

        self.time_col = time_col
        if data_cols:
            self.data_cols = data_cols
        else:
            self.data_cols = [i for i in range(self.n_cols) 
                              if i != self.time_col]

        # get dt and initial time from auxiliary data (if time included)
        if self.time_col is not None and self.constant_dt:
            self._initial_time = self.time
            self._read_next_step()
            self._dt = self.time - self._initial_time
            self.go_to_first_step()

    def __len__(self):
        return self.n_steps

    def next(self):
        """ Move to next step of auxiliary data """
        return self._read_next_step()

    def __next__(self):
        """ Move to next step of auxiliary data """
        return self.next()

    def __iter__(self):
        self._restart()
        return self

    def _restart(self):
        """ Reset back to start; calling next should read in first step """
        # Overwrite as appropriate
        self.step = -1
                
    def go_to_first_step(self):
        """ Return to and read first step """
        ## May need to overwrite, depending e.g. how header dealt with
        ## Could also use go_to_step(0) ...
        self._restart()
        self._read_next_step()
                
    def _read_next_step(self):  
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxilairy reader!")

    def read_ts(self, ts):
        """ Read auxiliary data for the *ts*.

        Read the auxiliary steps 'assigned' to *ts* (the steps that are within
        *ts.dt*/2 of of the trajectory timestep/frame - ie. closer to *ts*
        than either the preceeding or following frame). Then calculate a 
        'representative value' for the timestep from the data in each of these 
        auxiliary steps, and add to *ts*.

        Note
        ----
        The auxiliary reader will end up positioned at the last step assigned
        to the trajectory frame or, if the frame includes no auxiliary steps,
        (as when auxiliary data is less frequent), the most recent auxiliary 
        step before the frame.

        """
        # Make sure our auxiliary step starts at the right point (just before
        # the frame being read): the current step should be assigned to a 
        # previous frame, and the next step to either the frame being read of a 
        # following frame. Move to right position if not.
        if not (self.step_to_frame(self.step, ts) < ts.frame
                and self.step_to_frame(self.step+1, ts) >= ts.frame):
            return self.go_to_ts(ts)

        self.reset_ts() # clear previous ts data
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self.add_step_to_ts(ts.time)
        self.ts_rep = self.calc_representative()
        ts.aux.__dict__[self.name] = self.ts_rep
        return ts


    def step_to_frame(self, step, ts):
        """ Calculate closest trajectory frame for auxiliary step *step*.

        Calculated given dt and offset from *ts* as::
            math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)
        """
        if step >= self.n_steps or step < 0:
            ## make sure step is in the valid range. Raise error?
            return None 
        offset = ts.data.get('time_offset', 0)
        return math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)

    def go_to_ts(self, ts):
        """ Move to and read auxiliary steps corresponding to *ts* """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override go_to_ts() in auxiliary reader!")

    def reset_ts(self):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def add_step_to_ts(self, ts_time):
        """ Add data from the current step to *ts_data* """
        if len(self.ts_data) == 0:
            self.ts_data = np.array([self.step_data])
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts_time))

    def calc_representative(self):
        """ Calculate a represenatative auxiliary value(s).
        
        Currently available options for calculating represenatative value are
          *'closest': default; the value(s) from the step closest to in time to 
           the trajectory timestep
          *'average': average of the value(s) from steps 'assigned' to the 
           trajectory timestep. 
        Additionally, if a `cutoff` is specified, only steps within this time 
        of the trajectory timestep are considered in calculating the 
        represenatative.

        Returns
        -------
        List (of length n_cols) of auxiliary value(s) 'representing' the 
        timestep.
        """

        if self.cutoff != -1:
            cutoff_data = np.array([self.ts_data[i] 
                                    for i,d in enumerate(self.ts_diffs)
                                    if d < self.cutoff])
            cutoff_diffs = [d for d in self.ts_diffs if d < self.cutoff]
        else:
            cutoff_data = self.ts_data
            cutoff_diffs = self.ts_diffs
        if len(cutoff_data) == 0:
            value = np.array([np.nan]*len(self.data_cols))
        elif self.represent_ts_as == 'closest':
            value = cutoff_data[np.argmin(cutoff_diffs),:]
        elif self.represent_ts_as == 'average':
            value = np.mean(cutoff_data, axis=0)
        return value
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        # Overwrite as appropriate
        pass

    @property
    def time(self):
        """ Time in ps of current auxiliary step.

        As read from the appropriate column of the auxiliary data, if present; 
        otherwise calcuated as ``step * dt + initial_time``
        """
        if self.time_col is not None:
            return self._data[self.time_col]
        else:
            return self.step * self.dt + self.initial_time

    @property
    def step_data(self):
        """ Auxiliary values of interest for the current step.

        As taken from the appropariate columns (identified in `data_cols`) of 
        the full auxiliary data read in for the current step.
        """
        return [self._data[i] for i in self.data_cols]

    @property
    def n_steps(self):
        """ Total number of steps in the auxiliary data. """
        try:
            return self._n_steps
        except AttributeError:
            self._n_steps = self.count_n_steps()
            return self._n_steps

    @property
    def times(self):
        """ List of times of each step in the auxiliary data. 

        Calculated using `dt` and `initial_time` if `constant_dt` is True; 
        otherwise as read from each auxiliary step in turn. 
        """
          
        try:
            return self._times
        except AttributeError:
            if self.constant_dt:
                self._times = [i*self.dt+self.initial_time 
                               for i in range(self.n_steps)]
            else:
                self._times = self.read_all_times()
            return self._times

    @property
    def dt(self):
        """ Change in time between steps. 

        Defaults to 1ps if not provided or read from auxiliary data. """
        if self._dt:
            return self._dt
        else:
            return 1  ## default to 1ps; WARN?

    @property
    def initial_time(self):
        """ Time corresponding to first auxiliary step. 

        Defaults to 0ps if not provided or read from auxilairy data. """
        if self._initial_time:
            return self._initial_time
        else:
            return 0 ## default to 0; WARN?      


    
class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file.

    Extends AuxReader with methods particular to reading each step in turn
    from a file.
    """
    
    def __init__(self, auxname, filename, **kwargs):
        self.auxfilename = filename
        self.auxfile = open(filename)
        
        super(AuxFileReader, self).__init__(auxname, **kwargs)

    def close(self):
        """ Close auxfile, if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ Reposition to just before first step """
        self.auxfile.seek(0)
        self.step = -1
        
    def _reopen(self):
        self.auxfile.close()
        self.auxfile = open(self.auxfilename)
        self.step = -1
