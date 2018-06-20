
import unittest
import fusedwake.gcl.fortran as fgcl
import fusedwake.gcl.python.gcl as pygcl
import numpy as np
import os

from fusedwake.WindFarm import WindFarm
from fusedwake.gcl.interface import GCL


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
                               pygcl.get_r96(D, CT, TI, pars=[a1, a2, a3, a4, b1, b2]))

    def test_Rw(self):
        """Compare the two implementations of get_Rw
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
        R = D/2.
        x=D*np.linspace(0.,10.,100)

        self.assertTrue(np.array([np.allclose(fgcl.get_rw(x, D, CT, TI, a1, a2, a3, a4, b1, b2)[i],
                  pygcl.get_Rw(x, R, TI, CT, pars=[a1, a2, a3, a4, b1, b2])[i]) for i in range(3)]).all())

    def test_dU(self):
        """Compare the two implementations of get_dU
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
        R = D/2.
        x=D*np.linspace(0.,10.,100)
        r=D*np.linspace(0.,2.,100)

        self.assertTrue(np.allclose(fgcl.get_du(x,r,D,CT, TI, a1, a2, a3, a4, b1, b2),
                          pygcl.get_dU(x,r,R, CT, TI, pars=[a1, a2, a3, a4, b1, b2])))

    def test_dUeq(self):
        """Compare the two implementations of get_dUeq
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
        R = D/2.
        CT = 0.98
        TI = 0.10

        dx = np.array([6.,10.,15.])*D
        dy = np.array([1.,-5.,0.])*D
        dz = np.array([-2,2.,1.])*D

        # Wake operating turbines
        Rop = np.array([1.,2.,.5])*D
        Dop = 2.*Rop

        self.assertTrue(np.allclose(fgcl.get_dueq(dx,dy,dz,Dop,D,CT,TI,a1,a2,a3,a4,b1,b2),
                  pygcl.get_dUeq(dx,dy,dz,Rop,R,CT,TI,pars=[a1, a2, a3, a4, b1, b2])))

    def test_gcl(self):
        """Compare the two implementations of gcl
        """
        # Fixed parameters
        a1 = 0.435449861
        a2 = 0.797853685
        a3 = -0.124807893
        a4 = 0.136821858
        b1 = 15.6298
        b2 = 1.0

        current_dir = os.path.dirname(os.path.realpath(__file__))
        filename = 'test_WF_4Turbines.yml'
        yml_path = os.path.join(current_dir,filename)
        wf = WindFarm(yml=yml_path)
        gcl = GCL(WF=wf)

        # Inputs
        WD = np.arange(-50,50)+270.
        WS = 10.
        TI=0.1

        P_rat_py_v0 = []
        P_rat_py_v1 = []
        P_rat_fgcl = []
#        P_rat_fgclm_s = []
        for wd in WD:
            out = gcl(WF=wf, WS=WS, WD=wd, TI=TI, version='fort_gcl')
            P_rat_fgcl = np.append(P_rat_fgcl,out.p_wt[1]/out.p_wt[0])

#            out = gcl(WF=wf, WS=WS*np.ones([wf.nWT]),
#                  WD=wd*np.ones([wf.nWT]),
#                  TI=TI*np.ones([wf.nWT]), version='fort_gclm_s')
#            P_rat_fgclm_s = np.append(P_rat_fgclm_s,out.p_wt[1]/out.p_wt[0])

            out = gcl(WF=wf, WS=WS*np.ones([wf.nWT]),
                  WD=wd*np.ones([wf.nWT]),
                  TI=TI*np.ones([wf.nWT]),version='py_gcl_v1')
            P_rat_py_v1 = np.append(P_rat_py_v1,out.p_wt[1]/out.p_wt[0])

            out = gcl(WF=wf, WS=WS*np.ones([wf.nWT]),
                  WD=wd*np.ones([wf.nWT]),
                  TI=TI*np.ones([wf.nWT]), version='py_gcl_v0')
            P_rat_py_v0 = np.append(P_rat_py_v0,out.p_wt[1]/out.p_wt[0])

        WDm = WD.reshape([-1,1])*np.ones([1,wf.nWT])
        out = gcl(WF=wf, WS=WS*np.ones_like(WDm), WD=WDm,
                  TI=0.1*np.ones_like(WDm), version='fort_gclm')
        P_rat_fgclm = out.p_wt[:,1]/out.p_wt[:,0]

#        out = gcl(WF=wf, WS=WS*np.ones_like(WDm), WD=WDm,
#                  TI=0.1*np.ones_like(WDm), version='fort_gclm_av')
#        P_rat_fgclm_av = out.p_wt[:,1]/out.p_wt[:,0]

        self.assertTrue(np.allclose(P_rat_fgclm, P_rat_py_v0)&\
        np.allclose(P_rat_fgclm, P_rat_py_v1)&\
        np.allclose(P_rat_fgclm, P_rat_fgcl))#&\
#        np.allclose(P_rat_fgclm, P_rat_fgclm))

if __name__ == "__main__":
    unittest.main()
