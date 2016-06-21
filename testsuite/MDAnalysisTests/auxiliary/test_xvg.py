from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           raises)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import AUX_XVG, XVG_BAD_NCOL
from MDAnalysisTests.auxiliary.base import (BaseAuxReaderTest, BaseAuxReference)

class XVGReference(BaseAuxReference):
    def __init__(self):
        super(XVGReference, self).__init__()
        self.testdata = AUX_XVG
        self.reader = mda.auxiliary.XVG.XVGReader

class TestXVGReader(BaseAuxReaderTest):
    def __init__(self):
        reference = XVGReference()
        super(TestXVGReader, self).__init__(reference)

    def test_no_time_col(self):
        # XMGReader automatically sets time_col to 0 so we need to specifically
        # set it to none to test without
        self.reader = self.ref.reader(self.ref.testdata, dt=self.ref.dt, 
                                      initial_time=self.ref.initial_time,
                                      time_col=None)
        for i, val in enumerate(self.reader):
            assert_equal(val['data'], self.ref.all_data[i],
                         "step_data for step {0} does not match".format(i))

    @raises(ValueError)
    def test_wrong_n_col_raises_ValueError(self):
        # encountering a different number of columns at a later step should 
        # raise ValueError
        self.reader = self.ref.reader(XVG_BAD_NCOL)
        next(self.reader)
