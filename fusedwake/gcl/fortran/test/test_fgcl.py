
import unittest
from fusedwake.gcl.fortran import *
import numpy as np

current_dir = os.path.dirname(os.path.realpath(__file__))

class TestFortranGCL(unittest.TestCase):
    def setUp(self):
        pass

    def test_r96(self):
        # Fixed parameters
        a1 = 0.435449861
        a2 = 0.797853685
        a3 = -0.124807893
        a4 = 0.136821858
        b1 = 15.6298
        b2 = 1.0

        # Variables
        D = 80.0
        CT = 0.98
        TI = 0.10

        R96 = a1 * (np.exp(a2 * CT * CT+a3 * CT + a4)) * (b1 * TI + b2) * D
        self.assertAlmostEqual(get_r96(D, CT, TI, a1, a2, a3, a4, b1, b2), R96)

if __name__ == "__main__":
    unittest.main()
