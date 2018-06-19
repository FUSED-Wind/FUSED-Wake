#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_fusedwake
----------------------------------

Tests for `fusedwake` module.
"""

import unittest

from fusedwake import fusedwake


class TestFusedwake(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_000_something(self):
        pass

    def test_imports(self):
        import fusedwake.WindTurbine as wt
        import fusedwake.WindFarm as wf
        #import fusedwake.fused as fused_gcl
        from fusedwake.gcl.python.gcl import GCLarsen_v0


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
