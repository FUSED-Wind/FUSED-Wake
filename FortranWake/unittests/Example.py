import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import scipy as sp
import pandas as pd
from scipy.interpolate import pchip_interpolate as interp
from scipy.interpolate import pchip

import pickle

# Fused wrapped Fortran wake models
from FortranWake.fused_fortran import FNOJ, FGCL, AEP_f
from FortranWake.fusedwasp import PlantFromWWH, WTDescFromWTG

# FusedWind
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout, WTPC, WeibullWindRoseVT, GenericWindRoseVT

#------------------------------------

wt_layout = PlantFromWWH('hornsrev1_turbine_nodescription.wwh').wt_layout
P_rated = np.sum(wt_layout.wt_array('power_rating'))

#------------------------------------
#
# Improve the interpolation of the power curve and add CT_idle
#

for wt in wt_layout.wt_list:
    PCI = pchip(np.append(3.,wt.power_curve[:,0]),np.append(0.,wt.power_curve[:,1]))
    CTCI = pchip(np.append(3.,wt.c_t_curve[:,0]),np.append(0.053,wt.c_t_curve[:,1]))
    wt.cut_in_wind_speed = 3.
    u_c = np.linspace(3.,wt.cut_out_wind_speed,211,endpoint=True)
    p_c = PCI(u_c)
    ct_c = CTCI(u_c)
    wt.power_curve = np.vstack([u_c,p_c]).T
    wt.c_t_curve = np.vstack([u_c,ct_c]).T
    wt.c_t_idle = 0.053
    if wt.name[-2]=='_':
        wt.name='wt_0'+wt.name[-1]


df = pd.DataFrame([wt_layout.create_dict(n) for n in range(wt_layout.n_wt)])

#------------------------------------
#
# FNOJ
#

NOJ_s = FNOJ()
NOJ_s.wind_speeds=np.array(4.2)
NOJ_s.wind_directions= np.array(270.0)
NOJ_s.wake_expansions= np.array(0.05)
NOJ_s.wt_layout = wt_layout
NOJ_s.run()
NOJ_s.power/1e6

#------------------------------------
#
# FNOJ
#

wind_speeds = [i for i in np.linspace(4.0,25.,num=43,endpoint=True)]
wind_directions = [i for i in np.linspace(0.,360.,num=720,endpoint=False)]
WS,WD = np.meshgrid(np.array(wind_speeds),np.array(wind_directions))
WS,WD = WS.flatten(),WD.flatten()

NOJ = FNOJ()
NOJ.wind_speeds=WS
NOJ.wind_directions=WD
NOJ.wt_layout = wt_layout
NOJ.wake_expansions=0.050*np.zeros_like(WS)
NOJ.run()
NOJ.power/1e6 #Power[iWD,iWD]

#------------------------------------
#
# GCL
#

GCL_s = FGCL()
GCL_s.wind_speeds=np.array(4.2)
GCL_s.wind_directions=np.array(270.)
GCL_s.turbulence_intensities=np.array(0.07)
GCL_s.wt_layout = wt_layout
GCL_s.run()
GCL_s.power/1e6


#------------------------------------
#
# GCL
#

wind_speeds = [i for i in np.linspace(4.0,25.,num=22,endpoint=True)]
wind_directions = [i for i in np.linspace(0.,360.,num=91,endpoint=False)]
tis = [i for i in np.logspace(np.log10(0.04),np.log10(0.20),num=3,endpoint=True)]

WS,WD = np.meshgrid(np.array(wind_speeds),np.array(wind_directions))
WS,TI = np.meshgrid(WS.flatten(),np.array(tis))
WD,__ = np.meshgrid(WD.flatten(),np.array(tis))
WS,WD,TI = WS.flatten(),WD.flatten(),TI.flatten()

GCL = FGCL()
GCL.wind_speeds=WS
GCL.wind_directions=WD
GCL.turbulence_intensities=TI
GCL.wt_layout = wt_layout
GCL.run()
GCL.power/1e6 #Power[iWD,iWD]
