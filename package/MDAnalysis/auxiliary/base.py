import numpy as np
import math

class AuxReader(object):
    """ Base class for auxiliary readers
    Does iteration and aligning between trajectory ts/auxiliary step
    Actual reading/parsing of data to be handled by individual readers.
    """
      
    def __init__(self, auxnames, represent_ts_as='closest', cutoff=None, 
                 dt=None, initial_time=None, time_col=None, data_cols=None, 
                 constant_dt=True, **kwargs):

        self.names = auxnames
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff ## UNITS?
        # set initially to avoid error on first read
        self.n_cols = None

        self.ts_data = None
        self.ts_rep = None

        self._initial_time = initial_time
        self._dt = dt
        self.constant_dt=True

        self.time_col = time_col
        self.data_cols = data_cols

        self.step = -1
        self._read_next_step()
        self.n_cols = len(self._data)

        if self.time_col is not None and self.constant_dt:
            self._initial_time = self.time
            self._read_next_step()
            self._dt = self.time - self._initial_time
            self.go_to_first_step()

        if time_col >= self.n_cols:
            raise ValueError("Index {0} for time column out of range (num. "
                             "cols is {1})".format(time_col, self.n_cols))
        if data_cols:
            for col in data_cols:
                if col >= self.n_cols:
                    raise ValueError("Index {0} for data column out of range (num."
                                     " cols is {1})".format(col, self.n_cols))

        # update dt/initial time from data

    def next(self):
        """ Move to next step of data """
        return self._read_next_step()

    def __next__(self):
        """ Move to next step of data """
        return self.next()

    def __iter__(self):
        self._restart()
        return self

    def _restart(self):
        """ Reset back to start; calling next should read in first step """
        # Overwrite when reading from file to seek to start of file
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
        """ Read and record data from steps closest to *ts*, then 
        calculate representative value for ts """
        if not (self.step_to_frame(self.step, ts) < ts.frame
                and self.step_to_frame(self.step+1, ts) >= ts.frame):
            return self.go_to_ts(ts)

        self.reset_ts()
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self.add_step_to_ts(ts.time)
        self.ts_rep = self.calc_representative()
        ts.aux.__dict__[self.names] = self.ts_rep
        return ts
        ## after reading timestep, aux should be a last step in that timestep
        ## (or last step before that, if no steps in timestep)

    def step_to_frame(self, step, ts):
        if step >= self.n_steps or step < 0:
            return None ## raise error?
        offset = ts.data.get('time_offset', 0)
        return math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)

    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override go_to_ts() in auxilairy reader!")

    def reset_ts(self):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def add_step_to_ts(self, ts_time):
        if len(self.ts_data) == 0:
            self.ts_data = [self.step_data]
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts_time))

    def calc_representative(self):
        if self.cutoff:
            cutoff_data = np.array([self.ts_data[i] 
                                    for i,d in enumerate(self.ts_diffs)
                                    if d < self.cutoff])
            cutoff_diffs = [d for d in self.ts_diffs if d < self.cutoff]
        else:
            cutoff_data = self.ts_data
            cutoff_diffs = self.ts_diffs
        if len(cutoff_data) == 0:
            value = [] # TODO - flag as missing
        elif self.represent_ts_as == 'closest':
            value = cutoff_data[np.argmin(cutoff_diffs),:]
        elif self.represent_ts_as == 'average':
            value = np.mean(cutoff_data, axis=0)
        # TODO - interpolation?
        return value
    
    def close():
        # Overwrite when reading from file to close open file
        pass    

    @property
    def time(self):
        """ Time in ps of current auxiliary step
        As read from the auxiliary data, if present; otherwise calcuated
        as step * dt + initial_time
        """
        if self.time_col is not None:
            return self._data[self.time_col]
        else:
            return self.step * self.dt + self.initial_time
        ## warn/error if wrong number of columns?

    @property
    def step_data(self):
        if self.data_cols:
            return [self._data[i] for i in self.data_cols]
        else:
            return [self._data[i] for i in range(self.n_cols) 
                    if i != self.time_col]
        ## warn/error if wrong number of columns?

    @property
    def n_steps(self):
        try:
            return self._n_steps
        except AttributeError:
            self._n_steps = self.count_n_steps()
            return self._n_steps

    @property
    def times(self):
        if self.constant_dt:
            return [i*self.dt+self.initial_time for i in range(self.n_steps)]
        try:
            return self._times
        except AttributeError:
            self._times = self.read_all_times()
            return self._times

    @property
    def dt(self):
        if self._dt:
            return self._dt
        else:
            return 1  ## default to 1ps; WARN?

    @property
    def initial_time(self):
        if self._initial_time:
            return self._initial_time
        else:
            return 0 ## default to 0; WARN?      


    # TODO - add __enter__ and __exit__ methods when reading from file

    
class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file 
    Extends AuxReader with methods particular to reading from file"""
    
    def __init__(self, auxname, filename, **kwargs):
        self.auxfilename = filename
        self.auxfile = open(filename)
        
        super(AuxFileReader, self).__init__(auxname, **kwargs)

    def close(self):
        """ close file if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ reposition to just before first step """
        self.auxfile.seek(0)
        self.step = -1
        
    def _reopen(self):
        self.auxfile.close()
        self.auxfile = open(self.auxfilename)
        self.step = -1


