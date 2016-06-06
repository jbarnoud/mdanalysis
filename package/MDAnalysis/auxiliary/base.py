import numpy as np

class AuxReader(object):
    """ Base class for auxiliary readers 
    To read in from file or e.g. array
    """
      
    # TODO deal with changing/different units

    def __init__(self, auxdata, auxnames, represent_ts_as='closest', 
                 cutoff=None, dt=None, time_first_col=False, **kwargs):
        # TODO ? default auxnames; allow to name each column + pick which to keep
        # TODO pass in initial_time for when time not in data

        self.names = auxnames
        self.represent_ts_as = represent_ts_as

        self.step = 0
        
        self.cutoff = cutoff ## UNITS?

        self.ts_data = None
        self.ts_rep = None

        if dt:
            self.dt = dt
        elif time_first_col:
            # TODO set dt as diff between first two entries;
            # assumes dt is constant
        else:
            # TODO raise error - error must provide dt
        


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
        

