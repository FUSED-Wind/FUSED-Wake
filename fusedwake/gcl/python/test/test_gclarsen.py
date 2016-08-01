
import unittest
#from fusedwind.plant_flow.asym import AEPMultipleWindRoses
# from gclarsen.fused import FGCLarsen, FWindTurbine, FWindFarm, \
#     generate_GenericWindTurbinePowerCurveVT, generate_GenericWindFarmTurbineLayout, \
#     rosetta
import fusedwake.WindFarm as wf
import fusedwake.WindTurbine as wt
import fusedwake.gcl.python as gcl
import numpy as np
from fusedwake.py4we.wasp import WWH

import os

script_dir = os.path.dirname(os.path.realpath(__file__))

#from fusedwake.fusedwasp import PlantFromWWH
#
# class FGCLarsenTestCase(unittest.TestCase):
#
#     def setUp(self):
#         self.v80 = wt.WindTurbine('Vestas v80 2MW offshore','V80_2MW_offshore.dat',70,40)
#         self.HR1 = wf.WindFarm('Horns Rev 1','HR_coordinates.dat',self.v80)
#
#     def tearDown(self):
#         pass
#
#     def test_Init(self):
#         gcl = FGCLarsen()
#
#     def test_FWindTurbine(self):
#         """Generate a GenericWindTurbinePowerCurveVT from a WindTurbine, and
#         converte it back into a WindTurbine. Test that the conversion has not
#         corrupted any data.
#         """
#         fwt = FWindTurbine(self.v80.name, generate_GenericWindTurbinePowerCurveVT(self.v80))
#         # Test that the two turbines are identical
#         for i in ['H', 'R', 'u_cutin', 'u_cutout', 'name']:
#             self.assertEqual(getattr(fwt, i), getattr(self.v80, i))
#         for i in ['ref_P', 'ref_u', 'ref_CT']:
#             np.testing.assert_almost_equal(getattr(fwt, i), getattr(self.v80, i))
#
#     def test_FWindFarm(self):
#         fwf = FWindFarm(self.HR1.name, generate_GenericWindFarmTurbineLayout(self.HR1))
#         # Test that the two wind farms are identical
#         for i in ['name', 'nWT']:
#             self.assertEqual(getattr(fwf, i), getattr(self.HR1, i))
#         for i in ['name', 'nWT']:
#             self.assertEqual(getattr(fwf, i), getattr(self.HR1, i))
#
#
#     def test_Run(self):
#         inputs = dict(
#             WS=8.0,
#             z0=0.0001,
#             TI=0.05,
#             WD=270,
#             WF=self.HR1,
#             NG=4,
#             sup='lin',
#             pars=[0.5,0.9,-0.124807893,0.136821858,15.6298,1.0])
#         P_WT,U_WT, Ct = gcl.GCLarsen(**inputs)
#         fgcl = FGCLarsen()
#         # Setting the inputs
#         for k,v in rosetta.iteritems():
#             setattr(fgcl, v, inputs[k])
#         fgcl.wt_layout = generate_GenericWindFarmTurbineLayout(inputs['WF'])
#         fgcl.run()
#         np.testing.assert_almost_equal(P_WT, fgcl.wt_power)
#         np.testing.assert_almost_equal(U_WT, fgcl.wt_wind_speed)
#

class GCLarsen_v2_TestCase(unittest.TestCase):
    def setUp(self):
        self.v80 = wt.WindTurbine('Vestas v80 2MW offshore', script_dir+'/V80_2MW_offshore.dat', 70,40)
        self.HR1 = wf.WindFarm(name='Horns Rev 1', coordFile=script_dir+'/HR_coordinates.dat', WT=self.v80)
        self.inputs = dict(
            WS=8.0,
            z0=0.0001,
            TI=0.05,
            WD=270,
            WF=self.HR1,
            NG=4,
            sup='lin',
            pars=[0.5, 0.9, -0.124807893, 0.136821858, 15.6298, 1.0])

    def tearDown(self):
        pass

    def test_GCLarsen_v2(self):
       """Testing that the new implementation of GCLarsen is compatible with the
       old one.
       """#


       P_WT, U_WT, Ct = gcl.GCLarsen_v0(**self.inputs)
       P_WT2, U_WT2, Ct2 = gcl.GCLarsen(**self.inputs)

       np.testing.assert_almost_equal(P_WT, P_WT2)
       np.testing.assert_almost_equal(U_WT, U_WT2)
       np.testing.assert_almost_equal(Ct, Ct2)

# class test_AEP(unittest.TestCase):
#     def test_HR(self):
#         ### Single wind rose type
#         hrAEP = AEPMultipleWindRoses()
#         hrAEP.add('wf', FGCLarsen())
#         hrAEP.configure()
#         hrAEP.connect('wt_layout', 'wf.wt_layout')
#         hrAEP.wt_layout = PlantFromWWH(filename='wind_farms/horns_rev/hornsrev1_turbine_nodescription.wwh').wt_layout
#         hrAEP.wind_directions = [0., 90., 180., 270.]#linspace(0.0, 360.0, 3)[:-1]
#         hrAEP.wind_speeds = [8., 12., 24.]#linspace(4.0, 25.0, 3)
#         hrAEP.run()
#         print hrAEP.net_aep
#         print hrAEP.wt_aep



if __name__ == "__main__":
    unittest.main()
