
import unittest
import fusedwake.gcl.fortran as fgcl
import fusedwake.gcl.python as gcl
import numpy as np
import os


current_dir = os.path.dirname(os.path.realpath(__file__))

class TestGCLImplementations(unittest.TestCase):
    def setUp(self):
        pass

    def test_r96(self):
        """Compare the two implementations of R96
        """
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


        self.assertAlmostEqual(fgcl.get_r96(D, CT, TI, a1, a2, a3, a4, b1, b2),
                                gcl.get_r96(D, CT, TI, pars=[a1, a2, a3, a4, b1, b2]))

if __name__ == "__main__":
    unittest.main()
