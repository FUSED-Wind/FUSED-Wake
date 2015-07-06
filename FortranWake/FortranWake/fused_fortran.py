from numpy.testing import assert_equal

__all__ = ['FNOJensen']

from fusedwind.plant_flow.comp import GenericWindFarm
from fusedwind.interface import implement_base
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout, WTPC

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, VarTree, Array, List, Int, Enum, Str

# from topfarm.tlib import TopfarmComponent

import NOJ as f_noj #fortran implementation of NOJ
import GCL as f_gcl #fortran implementation of GCL

import numpy as np
from numpy import zeros,ones,ones_like,meshgrid,array,min,hstack,trapz,exp
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

@implement_base(GenericWindFarm)
class FusedFGCL(Component):
    """
    Implementation of the Gunner C. Larsen stationary wake model in according
    to fusedwind.plant_flow interface.
    """

    # Interface from GenericWindFarm
    # Inputs:
    wind_speed = Float(iotype='in', low=0.0, high=100.0, units='m/s',
        desc='Rotor averaged Inflow wind speed')
    wind_direction = Float(iotype='in', low=0.0, high=360.0, units='deg',
        desc='Inflow wind direction at hub height')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    pars = List([0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0],
        iotype='in', desc='GCLarsen model parameters')
    turbulence_intensity = Float(0.05, iotype='in', low='0.0', high='1.0',
        desc='Ambient turbulence intensity')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')

    # Outputs:
    power = Float(iotype='out', units='kW',
        desc='Total wind farm power production')
    thrust = Float(iotype='out', units='N',
        desc='Total wind farm thrust')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        """
        Fortran implementation of gcl model wraped for fusedwind
        """

        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)

        # Run the wind case
        self.wt_power, self.wt_wind_speed, self.wt_thrust = f_gcl.gcl(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speed,
            wd = self.wind_direction,
            ti = self.turbulence_intensity,
            a1 = self.pars[0],
            a2 = self.pars[1],
            a3 = self.pars[2],
            a4 = self.pars[3],
            b1 = self.pars[4],
            b2 = self.pars[5],
            ng = self.n_quad_points,
            rho = min(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed') )
            #ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        self.power = self.wt_power.sum()
        self.thrust = self.wt_thrust.sum()


@implement_base(GenericWindFarm)
class FusedFNOJ(Component):
    """
    Implementation of the N. O. Jensen stationary wake model in according
    to fusedwind.plant_flow interface.
    """

    # Interface from GenericWindFarm
    # Inputs:
    wind_speed = Float(iotype='in', low=0.0, high=100.0, units='m/s',
        desc='Inflow wind speed at hub height')
    wind_direction = Float(iotype='in', low=0.0, high=360.0, units='deg',
        desc='Inflow wind direction at hub height')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansion = Float(0.050, iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope')

    # Outputs:
    power = Float(iotype='out', units='kW',
        desc='Total wind farm power production')
    thrust = Float(iotype='out', units='N',
        desc='Total wind farm thrust')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        """
        Fortran implementation of noj model wraped for fusedwind
        """

        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind case
        self.wt_power, self.wt_wind_speed, self.wt_thrust = f_noj.noj(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = self.wind_speed,
            wd = self.wind_direction,
            kj = self.wake_expansion,
            rho = min(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed') )
            #ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        self.power = self.wt_power.sum()
        self.thrust = self.wt_thrust.sum()


class MultipleFusedFGCL(Component):

    # Inputs
    wind_speeds = List([], iotype='in', units='m/s',
        desc='The different wind speeds to run [n]')
    wind_directions = List([], iotype='in', units='deg',
        desc='The different wind directions to run [n]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    pars = List([0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0],
        iotype='in', desc='GCLarsen model parameters')
    turbulence_intensities = List([], iotype='in',
        desc='Ambient turbulence intensity to run [n]')
    n_quad_points = Int(4, low=4, high=6, iotype='in',
        desc='Number of points in the Gauss integration')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='The power production at each inflow [nWD, nWS]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind farm thrust')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):

        WS = array(self.wind_speeds)
        WD = array(self.wind_directions)
        TI = array(self.turbulence_intensities)

        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind case
        self.wt_power, self.wt_wind_speed, self.wt_thrust = f_gcl.gcl(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = WS,
            wd = WD,
            ti = TI,#self.turbulence_intensity*ones_like(WS),
            a1 = self.pars[0]*ones_like(WS),
            a2 = self.pars[1]*ones_like(WS),
            a3 = self.pars[2]*ones_like(WS),
            a4 = self.pars[3]*ones_like(WS),
            b1 = self.pars[4]*ones_like(WS),
            b2 = self.pars[5]*ones_like(WS),
            ng = self.n_quad_points,
            rho = min(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed') )
            #ct_idle = wt_layout.wt_array(attr='c_t_idle') )

        power = self.wt_power.sum(axis=1)
        thrust = self.wt_thrust.sum(axis=1)

        self.power  = power
        self.thrust  = thrust

class MultipleFusedFNOJ(Component):

    # Inputs
    wind_speeds = List([], iotype='in', units='m/s',
        desc='The different wind speeds to run [nWS]')
    wind_directions = List([], iotype='in', units='deg',
        desc='The different wind directions to run [nWD]')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='Wind turbine properties and layout')

    # Specific Inputs:
    wake_expansion = Float(0.050, iotype='in', low='0.0', high='1.0',
        desc='Linear wake radius expansion slope')

    # Outputs
    power = Array([], iotype='out', units='kW',
        desc='The power production at each inflow [nWD, nWS]')
    thrust = Array(iotype='out', units='N',
        desc='Total wind farm thrust')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):

        WS = array(self.wind_speeds)
        WD = array(self.wind_directions)

        # get T2T distance in global coordinates
        wt_layout = self.wt_layout
        x_g,y_g,z_g=get_T2T_gl_coord(wt_layout)
        # Run the wind case
        self.wt_power, self.wt_wind_speed, self.wt_thrust = f_noj.noj(
            x_g = x_g,
            y_g = y_g,
            z_g = z_g,
            dt = wt_layout.wt_array(attr='rotor_diameter'),
            p_c = wt_layout.wt_array(attr='power_curve'),
            ct_c = wt_layout.wt_array(attr='c_t_curve'),
            ws = WS,
            wd = WD,
            kj = self.wake_expansion*ones_like(WS),
            rho = min(wt_layout.wt_array(attr='air_density')),
            ws_ci = wt_layout.wt_array(attr='cut_in_wind_speed'),
            ws_co = wt_layout.wt_array(attr='cut_out_wind_speed') )
            #ct_idle = wt_layout.wt_array(attr='c_t_idle') )

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


    def __init__(self, wt_layout, wind_rose, wf, **kwargs):
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
        super(AEP_f, self).__init__(**kwargs)

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
        self.wf.wind_speeds=self.wind_speeds
        self.wf.wind_directions=self.wind_directions

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