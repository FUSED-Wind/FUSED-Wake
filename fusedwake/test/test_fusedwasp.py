from unittest import TestCase
from fusedwake.fusedwasp import  PlantFromWWH

__author__ = 'pire'


class TestPlantFromWWH(TestCase):

    def test__init__(self):
        wwh = PlantFromWWH(filename = 'wind_farms/horns_rev/hornsrev1_turbine_nodescription.wwh')
        assert (wwh.wt_layout.wt_names) > 0


    def test_call(self):
        pass
    def test_execute(self):
        pass
