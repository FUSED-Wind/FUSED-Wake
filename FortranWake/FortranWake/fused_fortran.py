from numpy.testing import assert_equal

__all__ = ['FNOJensen']

from fusedwind.plant_flow.comp import GenericWindFarm
from fusedwind.interface import implement_base
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout, WTPC

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, VarTree, Array, List, Int, Enum, Str

import NOJ as f_noj #fortran implementation of NOJ
import GCL as f_gcl #fortran implementation of GCL
import Mod_NOJ as f_mod_noj
import GAU as f_gau

import numpy as np
from numpy import zeros,ones,ones_like,meshgrid,array,hstack,trapz,exp,mean
from scipy.interpolate import interp1d

def get_T2T_gl_coord(wt_layout):
    """
    Function to calculated the turbine to turbine distances in the global
    coordinate system.

    Parameters
    ----------
    wt_layout   GenericWindFarmTurbineLayout (FusedWind)

    Returns
    -------
    x_g     x component of distance between Tj and Ti := x_g[i,j]
    y_g     y component of distance between Tj and Ti := y_g[i,j]
    z_g     z component of distance between Tj and Ti := z_g[i,j]

    """
    pos=wt_layout.wt_array(attr='position')
    # Compute the turbine to turbine vector in global coordinates
    x_g = zeros([wt_layout.n_wt,wt_layout.n_wt])
    y_g = zeros([wt_layout.n_wt,wt_layout.n_wt])
    z_g = zeros([wt_layout.n_wt,wt_layout.n_wt])
    for i in range(wt_layout.n_wt):
        for j in range(wt_layout.n_wt):
            x_g[i,j] = pos[j,0]-pos[i,0]
            y_g[i,j] = pos[j,1]-pos[i,1]
            z_g[i,j] = wt_layout.wt_array(attr='hub_height')[j] \
                     - wt_layout.wt_array(attr='hub_height')[i]

    return x_g,y_g,z_g

# ----------------------------------------------------------------------
#                FUSED-wind wrapped NOJ wake model
# ----------------------------------------------------------------------

class FNOJ(Component):
    """
    Implementation of the N. O. Jensen stationary wake model in according
    to fusedwind.plant_flow interface.
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c noj(x,y,z,DT,P_c,CT_c,WS,WD,kj)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c kj (float): Wake (linear) expansion coefficient
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_noj.noj(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            kj = self.wake_expansions,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


class FNOJ_AV(Component):
    """
    Implementation of the N. O. Jensen stationary wake model in according
    to fusedwind.plant_flow interface.

    This model includes the option of having wt not available
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')
    wt_available = Array([], iotype='in',
        desc='Defines if each wind turbine is available [nF, nWT]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='The power production at each inflow [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind farm thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine [nF, nWT]')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine [nF, nWT]')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c noj_av(x,y,z,DT,P_c,CT_c,WS,WD,kj,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c kj (float): Wake (linear) expansion coefficient
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power,self.wt_thrust,self.wt_wind_speed=f_noj.noj_av(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            kj = self.wake_expansions,
            av = self.wt_available,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust

# ----------------------------------------------------------------------
#                FUSED-wind wrapped GCL wake model
# ----------------------------------------------------------------------

class FGCL(Component):
    """
    Implementation of the G. C. Larsen stationary wake model in according
    to fusedwind.plant_flow interface.
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    turbulence_intensities = Array([], iotype='in',
        desc='Ambient turbulence intensities to run [nF]')
    pars = List([0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0],
        iotype='in', desc='GCLarsen model parameters')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c gcl(x,y,z,DT,P_c,CT_c,WS,WD,TI)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Rethore. A simple stationary semi-analytical
c wake model. Technical Report Riso, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Ambient turbulence intensity [-]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_gcl.gcl(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            ti = self.turbulence_intensities,
            a1 = self.pars[0]*ones_like(self.wind_speeds),
            a2 = self.pars[1]*ones_like(self.wind_speeds),
            a3 = self.pars[2]*ones_like(self.wind_speeds),
            a4 = self.pars[3]*ones_like(self.wind_speeds),
            b1 = self.pars[4]*ones_like(self.wind_speeds),
            b2 = self.pars[5]*ones_like(self.wind_speeds),
            ng = self.n_quad_points,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


class FGCL_AV(Component):
    """
    Implementation of the G. C. Larsen stationary wake model in according
    to fusedwind.plant_flow interface.

    This model includes the option of having wt not available
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    pars = List([0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0],
        iotype='in', desc='GCLarsen model parameters')
    turbulence_intensities = Array([], iotype='in',
        desc='Ambient turbulence intensity to run [nF]')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')
    wt_available = Array([], iotype='in',
        desc='Defines if each wind turbine is available [nF, nWT]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c gcl_av(x,y,z,DT,P_c,CT_c,WS,WD,TI,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Rethore. A simple stationary semi-analytical
c wake model. Technical Report Riso, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Ambient turbulence intensity [-]
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_gcl.gcl_av(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            ti = self.turbulence_intensities,
            av = self.wt_available,
            a1 = self.pars[0]*ones_like(self.wind_speeds),
            a2 = self.pars[1]*ones_like(self.wind_speeds),
            a3 = self.pars[2]*ones_like(self.wind_speeds),
            a4 = self.pars[3]*ones_like(self.wind_speeds),
            b1 = self.pars[4]*ones_like(self.wind_speeds),
            b2 = self.pars[5]*ones_like(self.wind_speeds),
            ng = self.n_quad_points,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


# ----------------------------------------------------------------------
#                FUSED-wind wrapped Modified NOJ wake model
# ----------------------------------------------------------------------

class FMod_NOJ(Component):
    """
    Implementation of the N. O. Jensen stationary wake model in according
    to fusedwind.plant_flow interface.
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c mod_noj(x,y,z,DT,P_c,CT_c,WS,WD,kj)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c kj (float): Wake (linear) expansion coefficient
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_mod_noj.mod_noj(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            kj = self.wake_expansions,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


class FMod_NOJ_AV(Component):
    """
    Implementation of the N. O. Jensen stationary wake model in according
    to fusedwind.plant_flow interface.

    This model includes the option of having wt not available
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')
    wt_available = Array([], iotype='in',
        desc='Defines if each wind turbine is available [nF, nWT]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='The power production at each inflow [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind farm thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine [nF, nWT]')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine [nF, nWT]')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c mod_noj_av(x,y,z,DT,P_c,CT_c,WS,WD,kj,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c kj (float): Wake (linear) expansion coefficient
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power,self.wt_thrust,self.wt_wind_speed=f_mod_noj.mod_noj_av(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            kj = self.wake_expansions,
            av = self.wt_available,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


# ----------------------------------------------------------------------
#                FUSED-wind wrapped GAU wake model
# ----------------------------------------------------------------------

class FGAU(Component):
    """
    Implementation of the Gaussian deficit stationary wake model in according
    to fusedwind.plant_flow interface.
    [1] Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model
    for wind-turbine wakes. Renewable Energy, 70, 116-123.
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c gau(x,y,z,DT,P_c,CT_c,WS,WD,ks)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c ks (array): Linear wake expansion [-]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_gau.gau(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            ks = self.wake_expansions,
            ng = self.n_quad_points,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


class FGAU_AV(Component):
    """
    Implementation of the G. C. Larsen stationary wake model in according
    to fusedwind.plant_flow interface.

    This model includes the option of having wt not available
    """
    # Inputs
    wind_speeds = Array([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nF]')
    wind_directions = Array([], iotype='in', units='deg',
        desc='The different wind directions to run [nF]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansions = Array([], iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope [nF]')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')
    wt_available = Array([], iotype='in',
        desc='Defines if each wind turbine is available [nF, nWT]')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='Total wind plant power production [nF]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind plant thrust [nF]')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        '''
c ----------------------------------------------------------------------
c gau_av(x,y,z,DT,P_c,CT_c,WS,WD,ks,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c ks (array): Linear wake expansion [-]
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
        '''
        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind flow case
        self.wt_power, self.wt_thrust, self.wt_wind_speed = f_gau.gau_av(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speeds,
            wd = self.wind_directions,
            ks = self.wake_expansions,
            av = self.wt_available,
            ng = self.n_quad_points,
            rho = mean(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed'),
            ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust


class AEP_f(Component):

    # Inputs
    wind_speeds = List([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nWS]')

    wind_directions = List([], iotype='in', units='deg',
        desc='The different wind directions to run [nWD]')

    wt_positions = Array(iotype='in')

    scaling = Float(1.0, iotype='in', desc='Scaling of the AEP')

    # Outputs
    array_aep = Array([], iotype='out', units='kW*h',
        desc='The energy production per sector [nWD, nWS]')

    gross_aep = Float(iotype='out', units='kW*h',
        desc='Gross Annual Energy Production before availability and loss impacts')

    net_aep = Float(iotype='out', units='kW*h',
        desc='Net Annual Energy Production after availability and loss impacts')

    capacity_factor = Float(0.0, iotype='out',
        desc='Capacity factor for wind plant')


    def __init__(self, wind_speeds, wind_directions, wt_positions, scaling,
        wt_layout, wind_rose, wf,**kwargs):
        """
        :param wt_layout: GenericWindFarmTurbineLayout()
        :param wind_rose: WeibullWindRoseVT()
        :param wf: MultipleFused wake model
        :param scaling: float [default = 1.0]
                        The scaling used to calculate the net_aep. If it is set to 0.0, the scaling
                        will be set to the net_aep the first time the simulation is run.
        """
        self.wf = wf
        self.wf.wt_layout = wt_layout
        self.wind_rose = wind_rose

        WS,WD = np.meshgrid(np.array(self.wind_speeds),
                np.array(self.wind_directions))
        WS,WD = WS.flatten(),WD.flatten()

        self.wf.wind_speeds=WS
        self.wf.wind_directions=WD
        super(AEP_f,self).__init__(**kwargs)

    def execute(self):
        power_curve = interp1d(
            np.mean(self.wf.wt_layout.wt_array(attr='power_curve'),axis=0)[:,0],
            np.sum(self.wf.wt_layout.wt_array(attr='power_curve'),axis=0)[:,1])\
            (self.wind_speeds)

        # build the cdf vector of the wind speed for each wind rose wind direction sector
        cdfw = []
        for iwd, wd in enumerate(self.wind_rose.wind_directions):
            cdfw.append(weibullCDF(self.wind_speeds, self.wind_rose.A[iwd], self.wind_rose.k[iwd]))

        # calculate the probability in each wind direction sector, using the CDF of the wind rose wind direction
        cdfd0 = [sum(self.wind_rose.frequency[:i]) for i in range(len(self.wind_rose.frequency)+1)]
        wd = hstack([self.wind_rose.wind_directions, [360]])
        cdfd1 = interp1d(wd, cdfd0)(self.wind_directions)

        net_aep = 0.0
        gross_aep = 0.0
        cwd = 0
        net_aeps = zeros([len(self.wind_directions)])
        gross_aeps = zeros([len(self.wind_directions)])

        self.wf.wt_layout.wt_positions=self.wt_positions
        self.wf.run()

        for iwd, wd in enumerate(self.wind_directions):
            if cwd < len(self.wind_rose.wind_directions):
                while wd >= self.wind_rose.wind_directions[cwd+1] and cwd < len(self.wind_rose.wind_directions)-2:
                    # switching wind rose wind direction sector
                    cwd += 1
            powers = zeros([len(self.wind_speeds)])
            for iws, ws in enumerate(self.wind_speeds):
                powers[iws] = self.wf.power[iwd,iws]

            # Integrating over the wind speed CDF
            net_aeps[iwd] = trapz(powers, cdfw[cwd]) * 365.0 * 24.0
            gross_aeps[iwd] = trapz(power_curve, cdfw[cwd]) * 365.0 * 24.0 #* self.wt_positions.shape[0]

        # Integrating over the wind direction CDF
        net_aep = trapz(net_aeps, cdfd1)
        gross_aep = trapz(gross_aeps, cdfd1)

        self.capacity_factor = net_aep / gross_aep

        print self.scaling
        if self.scaling == 0.0:
            # The scaling has to be calculated
            self.scaling = net_aep

        self.net_aep = net_aep / self.scaling
        self.gross_aep = gross_aep

def weibullCDF(x,A,k):
    """
    Returns the CDF of a weibull distribution over a vector x

    :param x: list or ndarray[n]
    :param A: Weibull coefficient A
    :param k: Weibull coefficient k
    :return: CDF
    """
    return 1.0 - exp(-(array(x)/A)**k)

# Rosetta stone to convert pure fortran-python inputs into FUSED-Wind wrapper inputs
rosettaGCL = {
    'ws':'wind_speed',
    'wd':'wind_direction',
    'ti':'turbulence_intensity',
    'ng':'n_quad_points'}

rosettaNOJ = {
    'ws':'wind_speed',
    'wd':'wind_direction',
    'kj':'wake_expansion'}