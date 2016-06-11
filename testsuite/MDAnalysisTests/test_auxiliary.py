
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           assert_array_almost_equal)
import MDAnalysis as mda

from MDAnalysisTests.datafiles import AUX_XVG

class BaseAuxReference(object):
    def __init__(self):
        self.n_steps = 5
        self.n_cols = 2
        self.all_step_data = [[1], [2], [4], [8], [16]]
        self.all_times = [0, 1, 2, 3, 4]
        self.dt = 1
        self.initial_time = 0

        self.ts_lowf = mda.coordinates.base.Timestep(0, dt=2)
        self.ts_lowf.frame=1
        self.ts_lowf_data = [[2], [4]]
        self.ts_lowf_rep = [4] 
        self.ts_lowf_rep_average = [3]
        self.ts_lowf_last_step = 2

        self.ts_highf = mda.coordinates.base.Timestep(0, dt=0.5)
        self.ts_highf.frame=1
        self.ts_highf_data = []
        self.ts_highf_rep = []
        self.ts_highf_last_step = 0


class BaseAuxReaderTest(object):
    def __init__(self, reference):
        self.ref = reference
        self.reader = self.ref.reader('test', self.ref.testdata)

    def test_n_steps(self):
        assert_equal(self.reader.n_steps, self.ref.n_steps)

    def test_n_cols(self):
        assert_equal(self.reader.n_cols, self.ref.n_cols)

    def test_dt(self):
        assert_equal(self.reader.dt, self.ref.dt)

    def test_initial_time(self):
        assert_equal(self.reader.initial_time, self.ref.initial_time)

    def test_first_step(self):
        self.reader.go_to_first_step()
        assert_equal(self.reader.step_data, self.ref.all_step_data[0])
        assert_equal(self.reader.time, self.ref.all_times[0])

    def test_next_to_second_frame(self):
        reader = self.ref.reader('test', self.ref.testdata)
        reader.next()
        assert_equal(reader.step_data, self.ref.all_step_data[1])
        assert_equal(reader.time, self.ref.all_times[1])

    def test_go_to_step(self):
        ## (not all aux readers might have go_to_step)
        self.reader.go_to_step(3)
        assert_equal(self.reader.step_data, self.ref.all_step_data[3])
        assert_equal(self.reader.time, self.ref.all_times[3])

    def test_next_past_last_step(self):
        self.reader.go_to_step(self.reader.n_steps-1)
        assert_raises(StopIteration, self.reader.next)

    def test_iter(self):
        for i, val in enumerate(self.reader):
            ##temp - assume time first column; change __iter__ to return
            ## data only or add full ref_data line to ref?
            ref_data = [self.ref.all_times[i]]+self.ref.all_step_data[i]
            assert_equal(val, ref_data)

    def test_read_ts_lowf(self):
        ## split up the following?
        self.reader.read_ts(self.ref.ts_lowf)
        assert_equal(self.reader.ts_data, self.ref.ts_lowf_data)
        assert_almost_equal(self.reader.ts_rep, self.ref.ts_lowf_rep)
        assert_almost_equal(self.ref.ts_lowf.aux.test, self.ref.ts_lowf_rep)
        assert_equal(self.reader.step, self.ref.ts_lowf_last_step)

    def test_read_ts_highf(self):
        ## split?
        self.reader.read_ts(self.ref.ts_highf)
        assert_equal(self.reader.ts_data, self.ref.ts_highf_data)
        assert_almost_equal(self.reader.ts_rep, self.ref.ts_highf_rep)
        assert_almost_equal(self.ref.ts_highf.aux.test, self.ref.ts_highf_rep)
        assert_equal(self.reader.step, self.ref.ts_highf_last_step)

    def test_ref_as_average(self):
        reader = self.ref.reader('test', self.ref.testdata, represent_ts_as='average')
        reader.read_ts(self.ref.ts_lowf)
        assert_almost_equal(reader.ts_rep, self.ref.ts_lowf_rep_average)

    ##TODO - file specific tests - opening/closing?

class XVGReference(BaseAuxReference):
    def __init__(self):
        super(XVGReference, self).__init__()
        self.testdata = AUX_XVG
        self.reader = mda.auxiliary.xvg.XVGReader

class TestXVGReader(BaseAuxReaderTest):
    def __init__(self):
        reference = XVGReference()
        super(TestXVGReader, self).__init__(reference)
