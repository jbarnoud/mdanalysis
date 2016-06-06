import numpy as np

class AuxReader(object):
    """ Base class for auxiliary readers 
    To read in from file or e.g. array
    """
      
    # TODO deal with changing/different units

    def __init__(self, auxdata, auxnames, represent_ts_as='closest', 
                 cutoff=None, dt=None, time_first_col=False, initial_time=None):
        # TODO ? default auxnames; allow to name each column + pick which to keep

        self.names = auxnames
        self.represent_ts_as = represent_ts_as

        self.step = 0
        
        self.cutoff = cutoff ## UNITS?

        self.ts_data = None
        self.ts_rep = None

        self.read_next_step()

        # set initial time; default to 0
        if time_first_col:
            init_t = self.time
        elif initial_time:
            init_t = intial_time
        else:
            init_t = 0 # TODO - warn using default?
        self.initial_time = init_t

        # set dt (assuming constant!); default to 1 ps.
        if time_first_col:
            self.read_next_step()
            dt = self.time - self.initial_time
        elif dt:
            self.dt = dt
        else:
            init_t = 1 # TODO - warn using default; (check units!)
        self.dt = dt
        

    def __iter__(self):
        self.reopen() ## ...if array...?
        return self

    def next(self):
        """ Move to next step of data """
        return self.read_next_step()

    def __next__(self):
        """ Move to next step of data """
        return self.next()

    # def read_next_step():  
        # DEFINE FOR EACH

    # TODO - add __enter__ and __exit__ methods



class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers dealing with files """
    def __init__(self, filename, auxname, **kwargs):
        super(AuxFileReader, self).__init__(filename, auxname, **kwargs)

        self.auxfilename = filename
        self.auxfile = open(filename)

    def rewind(self):
        """ position back at start and read first step """
        self.reopen()
        self.read_next_step()

    def close(self):
        """ close file if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def reopen(self):
        """ reposition to just before first step """
        self.auxfile.close()  
        self.auxfile = open(self.auxfilename) 
        self.step = 0
        

