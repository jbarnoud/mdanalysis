import numpy as np

class XVGReader:
    def __init__(self, filename, auxname, method='closest'):
        # [units?]
        self.xvgfilename = filename
        self.xvgfile = open(filename)
        self.name = auxname # allow to set different names for each column?
                            # default value? aux-1, aux-2, etc...
        self.rep_method = method # naming...
        self.step = 0
        # assuming xvg has time in first column; no need for dt/offset...
      
        self.read_next_step() ## read to ts of trajectory aux is being added to?
        
    def __iter__(self):
        self.reopen()
        return self

    def next(self):
        return self.read_next_step()

    def __next__(self):
        return self.next()

    def read_next_step(self):
        line = self.xvgfile.readline()
        if line:
            while line[0] in ['#', '@']:
                line = self.xvgfile.readline()
            self.time = float(line.split()[0])
            self.step_data = [float(i) for i in line.split()[1:]]
            # check number of columns is as expected...
            self.step = self.step + 1
            return self.step_data
        else: # better way to know EOF...?
            self.reopen()
            raise StopIteration

    def go_to_step(self, i):
        self.reopen()
        while self.step != i:
            value = self.read_next_step()
        return value    
        # reset back to starting step?
        
    def reopen(self):
        self.xvgfile.close()  
        self.xvgfile = open(self.xvgfilename) 
        self.step = 0
        self.read_next_step()
        
    def read_next_ts(self, ts):
        # check current times match up! --> go_to_ts
        ts_data = []
        diffs = []
        while self.time - ts.time < ts.dt/2.: ## best way to do this?
            ts_data.append(self.step_data)
            diffs.append(abs(self.time - ts.time))
            self.read_next_step()
        self.ts_data = np.array(ts_data)
        if len(self.ts_data) == 0:
            value = [] #flag as missing
        elif self.rep_method == 'closest':
            value = self.ts_data[np.argmin(diffs),:] ## better way for this...
        elif self.rep_method == 'average':
            value = np.mean(self.ts_data, axis=0)
            
        #ts.aux.__dict__['self.name'] = value
        return value
