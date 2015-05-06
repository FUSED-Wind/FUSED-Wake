
import unittest

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import pchip_interpolate as interp
import pickle
import time

# import from OpenMDAO
# ---------------------
from openmdao.util.eggsaver import save

# import from Fusedwind
# ---------------------
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout,\
    GenericWindTurbinePowerCurveVT, WTPC, WeibullWindRoseVT, GenericWindRoseVT
from fusedwind.plant_flow.generate_fake_vt import generate_random_wt_layout
from fusedwind.plant_flow.asym import AEPMultipleWindRoses

# import from FortranWake (Fusedwind wrapped Fortran wake models)
# ---------------------
from FortranWake.fused_fortran import FusedFNOJ, MultipleFusedFNOJ, \
    FusedFGCL, MultipleFusedFGCL, AEP_f, get_T2T_gl_coord, \
    rosettaGCL, rosettaNOJ
from FortranWake.fusedwasp import PlantFromWWH, WTDescFromWTG

# import from FortranWake (Fortran wake models)
# ---------------------
from FortranWake import GCL as FortranGCL
from FortranWake import NOJ as FortranNOJ

class FGCLarsenTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Init(self):
        pass

    def test_Run(self):

        wt_layout = generate_random_wt_layout()
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        dt = wt_layout.wt_array(attr='rotor_diameter')

        # Interpolate power curves to have the same number of elements
        nu = 22
        p_c = np.array([[[np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i],
                interp(wt_layout.wt_array(attr='power_curve')[j][:,0],
                wt_layout.wt_array(attr='power_curve')[j][:,1],
                np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i])] \
                for i in range(nu)] for j in range(wt_layout.n_wt)])

        ct_c = np.array([[[np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i],
                interp(wt_layout.wt_array(attr='c_t_curve')[j][:,0],
                wt_layout.wt_array(attr='c_t_curve')[j][:,1],
                np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i])] \
                for i in range(nu)] for j in range(wt_layout.n_wt)])
        for iwt, wt in enumerate(wt_layout.wt_list):
            wt.power_curve = p_c[iwt,:,:]
            wt.c_t_curve = ct_c[iwt,:,:]

        rho = np.min(wt_layout.wt_array(attr='air_density'))
        ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed')
        ws_co = wt_layout.wt_array(attr='cut_out_wind_speed')
        a1,a2,a3,a4,b1,b2=[0.5,0.9,-0.124,0.13,15.63,1.0]
        pars=[a1,a2,a3,a4,b1,b2]
        ws=8.0
        wd=270.0
        ti=0.07
        ng=5
        inputs=dict(
            ws=ws,
            wd=wd,
            ti=ti,
            ng=ng
            )

        P_WT,U_WT, Ct = FortranGCL.gcl(
                        x_g,y_g,z_g,dt,p_c,ct_c,ws,wd,ti,
                        a1,a2,a3,a4,b1,b2,ng,rho,ws_ci,ws_co)

        fgcl = FusedFGCL()
        # Setting the inputs
        fgcl.wt_layout = wt_layout
        for k,v in rosettaGCL.iteritems():
            setattr(fgcl, v, inputs[k])
        fgcl.pars=pars
        fgcl.run()
        if  np.allclose(P_WT, fgcl.wt_power, rtol=1.e-5, atol=1e-7):
            save(wt_layout, 'failures/FGCLarsenTestCase_'+ \
                time.strftime('%d_%m_%Y__%H_%M')+'.p', \
                fmt=4, proto=-1, logger=None)

        np.testing.assert_allclose(P_WT, fgcl.wt_power, rtol=1.e-5, atol=1e-7)
        np.testing.assert_allclose(U_WT, fgcl.wt_wind_speed, rtol=1.e-5, atol=1e-7)

class NOJensenTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Init(self):
        pass

    def test_Run(self):

        wt_layout = generate_random_wt_layout()
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        dt = wt_layout.wt_array(attr='rotor_diameter')

        # Interpolate power curves to have the same number of elements
        nu = 22
        p_c = np.array([[[np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i],
                interp(wt_layout.wt_array(attr='power_curve')[j][:,0],
                wt_layout.wt_array(attr='power_curve')[j][:,1],
                np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i])] \
                for i in range(nu)] for j in range(wt_layout.n_wt)])

        ct_c = np.array([[[np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i],
                interp(wt_layout.wt_array(attr='c_t_curve')[j][:,0],
                wt_layout.wt_array(attr='c_t_curve')[j][:,1],
                np.linspace(wt_layout.wt_array(attr='cut_in_wind_speed')[j],
                wt_layout.wt_array(attr='cut_out_wind_speed')[j],nu)[i])] \
                for i in range(nu)] for j in range(wt_layout.n_wt)])
        for iwt, wt in enumerate(wt_layout.wt_list):
            wt.power_curve = p_c[iwt,:,:]
            wt.c_t_curve = ct_c[iwt,:,:]

        rho = np.min(wt_layout.wt_array(attr='air_density'))
        ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed')
        ws_co = wt_layout.wt_array(attr='cut_out_wind_speed')
        ws=8.0
        wd=270.0
        kj=0.05
        inputs=dict(
            ws=ws,
            wd=wd,
            kj=kj
            )

        P_WT,U_WT, Ct = FortranNOJ.noj(
                        x_g,y_g,z_g,dt,p_c,ct_c,ws,wd,kj,rho,ws_ci,ws_co)

        fnoj = FusedFNOJ()
        # Setting the inputs
        fnoj.wt_layout = wt_layout
        for k,v in rosettaNOJ.iteritems():
            setattr(fnoj, v, inputs[k])
        fnoj.run()
        np.testing.assert_almost_equal(P_WT, fnoj.wt_power)
        np.testing.assert_almost_equal(U_WT, fnoj.wt_wind_speed)




class test_AEP(unittest.TestCase):
    def test_HR1(self):
        ### Single wind rose type
        wwh=PlantFromWWH(filename='hornsrev1_turbine_nodescription.wwh')
        hr1_aep = AEP_f(wt_layout=wwh.wt_layout,
          wind_rose=getattr(wwh.wind_rose_vt,wwh.wt_layout.wt_names[0]),
          wf=MultipleFusedFGCL(),
          wind_speeds= np.linspace(4.0, 25.0, 10),#[8., 12., 24.],#
          wind_directions=np.linspace(0.0, 360.0, 15)[:-1],#[0., 90., 180., 270.],#
          scaling=1.0,
          wt_positions=wwh.wt_layout.wt_array(attr='position'),
          turbulence_intensity=0.07)

        hr1_aep.run()
        print
        print 'test_HR1'
        print hr1_aep.net_aep
        print hr1_aep.gross_aep
        print hr1_aep.capacity_factor



if __name__ == "__main__":
    unittest.main()
