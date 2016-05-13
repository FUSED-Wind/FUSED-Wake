from fusedwake.WindFarm import WindFarm
import unittest
import os
import numpy as np
current_dir = os.path.dirname(os.path.realpath(__file__))

class TestWindFarm(unittest.TestCase):
    def setUp(self):
        filename = current_dir + '/../../examples/middelgrunden.yml'
        self.wf = WindFarm(name='farm_name', yml=filename)

    def test_get_T2T_gl_coord(self):
        np.testing.assert_array_almost_equal(
            np.array(self.wf.get_T2T_gl_coord()),
            np.array(self.wf.get_T2T_gl_coord2()))

if __name__ == '__main__':
    unittest.main()
