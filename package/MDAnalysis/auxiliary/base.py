import numpy as np

class AuxReader(object):
    """ Base class for auxiliary readers
    Does iteration and aligning between trajectory ts/auxiliary step
    Actual reading/parsing of data to be handled by individual readers.
    """
      
    # TODO deal with changing/different units

    def __init__(self, auxnames, represent_ts_as='closest', cutoff=None, 
                 dt=None, initial_time=None, time_col=None, data_cols=None, **kwargs):

        self.names = auxnames
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff ## UNITS?

        self.ts_data = None
        self.ts_rep = None

        self._initial_time = initial_time
        self._dt = dt

        self.time_col = time_col
        self.data_cols = data_cols

        self.n_cols
        if time_col >= self.n_cols:
            raise ValueError("Index {0} for time column out of range (num. "
                             "cols is {1})".format(time_col, self.n_cols))
        if data_cols:
            for col in data_cols:
                if col >= self.n_cols:
                    raise ValueError("Index {0} for data column out of range (num."
                                     " cols is {1})".format(col, self.n_cols))

        self.go_to_first_step()

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
        self.step = 0
                
    def go_to_first_step(self):
        """ Return to and read first step """
        ## May need to overwrite, depending e.g. how header dealt with
        ## Could also use go_to_step(1) ...
        self._restart()
        self._read_next_step()
                
    def _read_next_step(self):  
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxilairy reader!")

    def read_ts(self, ts):
        """ Read and record data from steps closest to *ts*, then 
        calculate representative value for ts """
        # Make sure auxiliary and trajectory are still aligned
        #if ts.frame != 0 and not self.first_in_ts(ts):
         #   return self.go_to_ts(ts)
             ## TODO - this doesn't work when no aux for that step!!

        self.reset_ts()
        while self.step_in_ts(ts):
            self.add_step_to_ts(ts.time)
            try:
                self._read_next_step()
            except StopIteration:
                break
        self.ts_rep = self.calc_representative()
        ts.aux.__dict__[self.names] = self.ts_rep
        return ts
        ## currently means that after reading in a timestep, ts_data and
        ## ts_diffs correspond to that ts but the current step/step_data of 
        ## the auxreader is the first step 'belonging' of the next ts...

    def first_in_ts(self, ts):
        """ Check if current step is first step 'belonging' to *ts*
        Assumes auxiliary *dt* is constant! """
        if (self.time-(ts.time-ts.dt/2)) < self.dt:
            return True
        else:
            return False

    def step_in_ts(self, ts):
        """ Check if current step 'belongs' to *ts* """
        if (self.time-ts.time) <= ts.dt/2. and (self.time-ts.time) > -ts.dt/2.:
            return True
        else:
            return False
           
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
        as (step - 1) * dt + initial_time
        """
        if self.time_col is not None:
            return self._data[self.time_col]
        else:
            return (self.step - 1) * self.dt + self.initial_time
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
    def n_frames(self):
        try:
            return self._n_frames
        except AttributeError:
            self.get_info_from_data()
            return self._n_frames

    @property
    def n_cols(self):
        try:
            return self._n_cols
        except AttributeError:
            self.get_info_from_data()
            return self._n_cols

    @property
    def dt(self):
        if self._dt:
            return self._dt
        elif self.time_col is not None:
            self.get_info_from_data()
            return self._dt
        else:
            return 1  ## default to 1ps; WARN?

    @property
    def initial_time(self):
        if self._initial_time:
            return self._initial_time
        elif self.time_col is not None:
            self.get_info_from_data()
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
        self.step = 0
        
    def _reopen(self):
        self.auxfile.close()
        self.auxfile = open(self.auxfilename)
        self.step = 0
