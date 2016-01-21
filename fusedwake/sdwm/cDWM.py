# -*- coding: utf-8 -*-
""" Multiple classes definition from the DWM flowfield model
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
from math import pi

class Meta:
    """ Main class holding all meta data for sDWM core model
    """

    def __init__(self):
        """ Initializing the data / Set default value ."""

        ## Default Ambient conditions
        self.WS                    = 8.0 # wind speed m/s
        self.TI                    = 0.07 # turbulence intensity % X10 ?

        ## Default Park settings:
        self.WTG                   = 'V80'

        ## Domain settings:
        self.lx                    = 3.2 # global domain width in R
        self.ly                    = 3.2 # global domain height in R
        self.hub_x                 = self.lx/2. # position of the rotor center
        self.hub_y                 = self.ly/2.
        self.hub_z                 = np.asarray([]) # Turbine spacing [D]
        self.dR                    = 20 # points per R in radial direction default 20
        self.dx                    = 20 # points per R in lateral direction  default 20
        self.dy                    = 20 # points per R in longitudinal direction  default 20
        self.dz                    = 10 # points per D in axial direction default 10
        self.nx                    = self.lx * self.dx  # nb points in x (lateral)  direction global flow field
        self.ny                    = self.ly * self.dy  # nb points in x (lateral)  direction global flow field
        self.nz                    = np.asarray([])  # nb points in z (streamwise) direction global flow field
        self.x_vec                 = np.arange(1,self.nx+1,1)*(1.0/self.dx) # grid coordinates in lat. direction [R]
        self.y_vec                 = np.arange(1,self.ny+1,1)*(1.0/self.dy) # grid coordinates in long. direction [R]
        self.z_vec                 = np.asarray([]) # grid coordinates in streamwise. direction [D]
        self.z_vec_old             = np.asarray([]) # grid coordinates in streamwise. direction [D] at previous turbine
        self.x_mat                 = np.tile(self.x_vec.reshape(len(self.x_vec),1),self.ny)  # matrix grid holding x coordinates
        self.y_mat                 = np.tile(self.y_vec,(self.nx,1)) # matrix grid holding y coordinates
        self.z_mat                 = np.asarray([]) # matrix grid holding z coordinates
        self.nr                    = round(np.sqrt((max(self.x_vec - self.hub_x))**2 + (max(self.y_vec - self.hub_y))**2 ) * self.dR) # number of position in radial coordinate system
        self.r_dist_2              = np.sqrt((self.x_mat - self.hub_x)**2 + (self.y_mat - self.hub_y)**2 )  #Find distance to centre of wake plane
        self.C2C                   = np.asarray([]) #the "Center to Center" distance between hubs

        # Eddy viscosity calc mixL model specific parameter
        self.R_WTG                 = 1.0 # normalized turbine radius for eddy-viscosity calculation
        self.D_WTG                 = 2.0 * self.R_WTG # normalized turbine radius for eddy-viscosity calculation
        self.U0                    = 1.0 # normalized free stream velocity for eddy-viscosiy calculation
        self.dz_mixl               = self.D_WTG / self.dz  #note  is per D for eddy-viscosity calculation
        self.dr_mixl               = self.dR # points per R in radial direction identical to global flow field
        self.dr                    = self.R_WTG / self.dr_mixl  #note is per R for eddy-viscosity calculation
        self.lr_mixl               = self.lx+1. # domain size in radial direction [R] (Should be about 0.5R wider flow field domain width)
        self.vr_mixl               = np.linspace(0,self.lr_mixl-self.dr,(self.lr_mixl)/self.R_WTG*self.dr_mixl)  # coordinate vector
        self.vr_m                  = np.arange(1,self.dr_mixl*self.lr_mixl+1,1,dtype=float) / self.dr_mixl  # coordinate vector for polar discretization
        self.dA_DWM                = np.concatenate(([0], (pi*self.vr_m[1:]**2 - pi*self.vr_m[0:-1]**2)), axis=0)  # area vector for polar discretization
        self.lz_mixl               = np.asarray([]) # in R: 1R longer than global flow field due to backward diff scheme, turbine position dependent

        ## Model inputs calibrated constant
        self.Model                 = 'mixL' # Eddy viscosity model as per Keck et al. [2]
        # self.fU                    = 1.10 # fU = 1.1 & fR = 0.98 yields the best fit to AL calibration data according to Keck et al. [5].
        # self.fR                    = 0.98 #
        self.fU                    = 1. #Madsen et al [2].
        self.fR                    = 1.
        self.atmo_stab             = 'N'
        self.k1                    = 0.0919 # k1 = 0.0919 & k2 = 0.0178 yields the best fit to AL calibration data according to Keck et al. [5].
        self.k2                    = 0.0178

        ## Model flags
        self.Tbuildup_setting      = 1 # 1 = on # flag to control the buildup of turbulence along turbine rows, if disabled TI is equal to free stream TI
        self.wake_ind_setting      = 1 # 1 = on # flag to control the wake accumulation along turbine rows, if disabled U is equal to free stream velocity
        self.accu                  = 'dominant' # accumulation model
        self.accu_inlet            = True # accumulation of inlet velocity field. If enabled, the inflow flow field to each turbine is computed as the aggregated upstream flow field. If disabled, the model behaves like the HAWC2-DWM model [2] and [4]
        self.full_output           = False # control the amount of outputs scalar. If true, the velocity flow field is returned

        ## Misc variables
        self.mek_elec_power        = 0.95 # mek power to electric power
        self.rho                   = 1.25 # air at +10C
        self.wtg_ind               = np.asarray([]) # turbine index in main loop
        self.WTG_spec              = [] # class holding the turbine spec

        ## DWM Filter functions
        # Based on calibration data according to Keck et al. [5].
        self.f1                    =[0., 4.,] # F1 function = 1 after 4R
        self.f2                    =[0.035, 0.35] # F2 function starts at 0.035 for first 4R, then follows F2 = 1-0.965*e(-0.45*(2X-4)) (X in [R])

    def parse(self,**par):
        """Parsing data to class holding all meta data for sDWM core model."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')


        # here we update the attributes that are affected by the new variable parsing
        self.mean_WS=[]
        self.mean_TI=[]
        self.mean_WS_DWM=[]
        self.mean_TI_DWM=[]
        self.power_estimate=[]
        self.bem_power_estimate=[]

        # self.hub_x                 = self.lx/2. # position of the rotor center
        self.hub_y                 = self.ly/2.
        self.nx                    = self.lx * self.dx  # nb points in x (lateral)  direction global flow field
        self.ny                    = self.ly * self.dy  # nb points in x (lateral)  direction global flow field
        self.nz                    = np.asarray([])  # nb points in z (streamwise) direction global flow field
        self.x_vec                 = np.arange(1,self.nx+1,1)*(1.0/self.dx) # grid coordinates in lat. direction [R]
        self.y_vec                 = np.arange(1,self.ny+1,1)*(1.0/self.dy) # grid coordinates in long. direction [R]
        self.z_vec                 = np.asarray([]) # grid coordinates in streamwise. direction [D]
        self.z_vec_old             = np.asarray([]) # grid coordinates in streamwise. direction [D] at previous turbine
        self.x_mat                 = np.tile(self.x_vec.reshape(len(self.x_vec),1),self.ny)  # matrix grid holding x coordinates
        self.y_mat                 = np.tile(self.y_vec,(self.nx,1)) # matrix grid holding y coordinates
        self.z_mat                 = np.asarray([]) # matrix grid holding z coordinates
        self.lr_mixl               = self.lx+1. # 10 for single wake model, 6 for field model ? domain size in radial direction [R] (Should be about 0.5R wider flow field domain width)
        self.vr_mixl               = np.linspace(0,self.lr_mixl-self.dr,(self.lr_mixl)/self.R_WTG*self.dr_mixl)  # coordinate vector
        self.vr_m                  = np.arange(1,self.dr_mixl*self.lr_mixl+1,1,dtype=float) / self.dr_mixl  # coordinate vector for polar discretization
        self.dA_DWM                = np.concatenate(([0], (pi*self.vr_m[1:]**2 - pi*self.vr_m[0:-1]**2)), axis=0)  # area vector for polar discretization
        # self.nr                    = round(np.sqrt((max(self.x_vec - self.hub_x))**2 + (max(self.y_vec - self.hub_y))**2 ) * self.dR)
        # self.r_dist_2              = np.sqrt((self.x_mat - self.hub_x)**2 + (self.y_mat - self.hub_y)**2 )  #Find distance to centre of wake plane

class Meand:
    """ Class holding meandering magnitude related data
    """
    def __init__(self):
        """ Initializing the data / Set default value ."""

        ## Default Ambient conditions
        self.time                  = np.arange(1,600+1,1) # time vector in second for the meandering process
        self.std_meand_x           = np.asarray([]) # std deviation lateral of wake center
        self.std_meand_y           = np.asarray([]) # std deviation long of wake center
        self.meand_pos_x           = np.asarray([]) # wake position vector lat
        self.meand_pos_y           = np.asarray([]) # wake position vector long

    def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')

class FFoR:
     def __init__(self):
        """ Initializing the data / Set default value ."""
        self.x_vec_t               = np.asarray([])    #coordinate vector in fixed frame of reference nx,nz
        self.x_mat_t               = np.asarray([])    # coordinate matrix in fixed frame of reference (nx,ny,nz))

     def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')



class MFoR:
    """ Initialize structure holding calc results """

    def __init__(self,name):
        """ Initializing the data ."""
        self.name=name
        self.U=np.array([]) # streamwise velocity matrix in meandering frame of reference (vz_mixl,vr_mixl)
        self.V=np.array([] )# lateral velocity matrix in meandering frame of reference (vz_mixl,vr_mixl)
        self.visc=[]        # total viscosity matrix (vz_mixl,vr_mixl)
        self.TI_DWM=[]      # turbulence intensity forulated based on the relation derived between u'v' and u'u' (see Keck et al. [5])
        self.Turb_Stress_DWM=[] # matrix holding turbulence stress
        self.du_dr_tot=[]    # total velocity gradient (meandering + ABL)
        self.du_dr_DWM=[]    # velocity gradient from meandering
        self.Shear_add_du_dr=[] # added shear velocity gradient
        self.U_init =[]   # initial radial velocity deficit (axisymetric) input to ainslie model


    def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')


class Outputs:
    """ Initialize structure holding calc results """

    def __init__(self):
        """ Initializing the data ."""
        self.WS_axial_sym           = np.asarray([])
        self.WS_axial_sym_pow3      = np.asarray([])
        self.TI_axial_sym           = np.asarray([])
        self.TI_meand_axial_sym     = np.asarray([])
        self.TI_tot_axial_sym       = np.asarray([])
        # Store Mean turbine data
        self.C2C                     = np.asarray([])
        self.WS                     = np.asarray([])
        self.TI                     = np.asarray([])
        self.WS_DWM_BC              = np.asarray([])
        self.TI_DWM_BC              = np.asarray([])

        self.U_init_raw=np.array([])
        self.U_init=np.array([])
        self.vr_mixl=np.array([])
        self.vz_mixl=np.array([])

        self.dA=np.array([])
        self.mean_a=np.array([])
        self.f_w=np.array([])
        self.r_w=np.array([])
        self.U_w=np.array([])

        self.Shear_add_du_dr=np.array([])
        self.Shear_add_du_dz=np.array([])
        self.TI_DWM=np.array([])
        self.Turb_Stress_DWM=np.array([])
        self.U=np.array([])
        self.U_init=np.array([])
        self.V=np.array([])
        self.WakeW=np.array([])
        self.du_dr_DWM=np.array([])
        self.du_dr_tot=np.array([])
        self.visc=np.array([])


    def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')


class Aero:

    """ Load inductions from available libraries based on turbine codename """
    def __init__(self,name):
        """ Initializing the data ."""
        self.name=name
        self.a=np.array([])
        self.r_t=np.array([])
        self.CP=np.array([])
        self.CPloc=np.array([]) # local CP on BEM grid
        self.Power=np.array([])
        self.pow_cur=np.array([])
        self.ct_cur=np.array([])
        self.dA=[]
        self.mean_a=[]
        self.f_w=[]
        self.r_w=[]
        self.U_w=[]
        self.Un=[]
        self.Ut=[]
        self.CT=np.array([])
        self.Cp=[]  # local Cp on flow field grid at rotor
        self.RPM=[]
        self.PITCH=[]
        self.RPM_opt=[]
        self.PITCH_opt=[]