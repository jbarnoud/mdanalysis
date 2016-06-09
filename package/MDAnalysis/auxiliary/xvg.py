import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Read data from .xvg file
    
    Assumes data is time-ordered and first column is time

    Currently reading from file; can probably load a once and read through array
    instead
    """
 
    def __init__(self, auxname, filename, **kwargs):
        super(XVGReader, self).__init__(auxname, time_col=0, **kwargs)

        
    def _read_next_step(self):
        """ Read next recorded timepoint in xvg file """
        line = self.auxfile.readline()
        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self._data = [float(i) for i in line_no_comment.split()]
            if len(self._data) != self.n_cols:
                pass ## TODO - error?
            self.step = self.step + 1
            return self.step_data
        else:
            self.go_to_first_step()
            raise StopIteration
 
    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        self.go_to_first_step()
        while not self.step_in_ts(ts):
            self._read_next_step()
        return self.read_next_ts(ts)

    def go_to_step(self, i):
        """ Move to and read i-th step """
        ## probably not the best way to do this - seek?
        if not isintance(i, int):
            raise TypeError("Step number must be integer")
        if i > self.n_steps:
            raise ValueError("{0} is out of range range of auxiliary"
                             "(num. steps {1}!").format(i, self.n_steps))
        if i < 1:
            raise ValueError("Step numbering begins at 1")

        self.go_to_first_step()
        while self.step != i:
            value = self._read_next_step()
        return value

    def get_info_from_data(self):
        self._restart()
        self._read_next_step()
        if self.time_col is not None:
            self._initial_time = self.time
        self._n_cols = len(self._data)
        self._read_next_step()
        if self.time_col is not None:
            self._dt = self.time - self._initial_time
        i = 2
        while True:
            try:
                self._read_next_step()
                i = i + 1
            except StopIteration:
                break
        self._n_steps = i

