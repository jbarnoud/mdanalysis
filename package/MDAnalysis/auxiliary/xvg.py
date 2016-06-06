import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Read data from .xvg file
    
    Assumes data is time-ordered and first column is time
    """

    def __init__(self, filename, auxname, **kwargs):
        super(XVGReader, self).__init__(filename, auxname, **kwargs)

        self.read_next_step()
        
    def read_next_step(self):
        """ Read next recorded timepoint in xvg file """
        line = self.auxfile.readline()

        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self.time = float(line_no_comment.split()[0])
            self.step_data = [float(i) for i in line_no_comment.split()[1:]]
            # TODO check number of columns is as expected...
            self.step = self.step + 1
            return self.step_data
        else:
            self.reopen()
            raise StopIteration



    ## [the following can probably largely migrate to somewhere in base]
    
    def read_next_ts(self, ts):
        """ Read and record data from steps closest to *ts*
        Calculate representative value for ts """
        # Make sure auxiliary and trajectory are still aligned
        if not self.step_in_ts(ts) or not self.first_in_ts(ts):
            return self.go_to_ts(ts)
        else:
            self.reset_ts()
            while self.step_in_ts(ts):
                self.add_step_to_ts(ts)
                self.read_next_step()
            ## TODO - catch StopIteration error reach end of auxfile
            self.ts_rep = self.calc_representative()
            ts.aux.__dict__[self.names] = self.ts_rep
            return ts
        ## currently means that after reading in a timestep, ts_data and
        ## ts_diffs correspond to that ts but the current step/step_data of 
        ## the auxreader is the first step 'belonging' of the next ts...

    def reset_ts(self):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def step_in_ts(self, ts):
        """ Check if current step 'belongs' to *ts* """
        if (self.time-ts.time) <= ts.dt/2. and (self.time-ts.time) > -ts.dt/2.:
            return True
        else:
            return False

    def first_in_ts(self, ts):
        """ Check if current step is first step belonging to *ts*
        Assumes auxiliary *dt* is constant """
        if (self.time-ts.time+ts.dt/2) < self.dt/2:
            return True
        else:
            return False

    def add_step_to_ts(self, ts):
        if len(self.ts_data) == 0:
            self.ts_data = [self.step_data]
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts.time))
       
    def calc_representative(self):
        if self.cutoff:
            ## starting to get nasty, should probably change...
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
        return value

    def go_to_step(self, i):
        """ Move to and read i-th step """
        ## probably not the best way to do this - seek?
        self.reopen()
        while self.step != i:
            value = self.read_next_step()
        return value
        
    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        self.rewind()
        while not self.step_in_ts(ts):
            self.read_next_step()
        return self.read_next_ts(ts)
