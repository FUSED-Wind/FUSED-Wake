from numpy.testing import assert_equal

__all__ = ['FGCLarsen']

from fusedwind.plant_flow.comp import GenericWindFarm
from fusedwind.interface import implement_base
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout,\
    GenericWindTurbinePowerCurveVT, WTPC

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, VarTree, Array, List, Int, Enum, Str
from WindFarm import WindFarm
from WindTurbine import WindTurbine
import GCLarsen as gcl
import numpy as np

# Make sure that your class has some kind of docstring. Otherwise
# the descriptions for your variables won't show up in the
# source ducumentation.

class FWindTurbine(WindTurbine):
    def __init__(self, name, wt_desc, CT_idle=0.053):
        """Initializes a WindTurbine instance from a GenericWindTurbinePowerCurveVT
        instance

        Inputs:
           name (str):  Wind turbine name
           wt_desc (GenericWindTurbinePowerCurveVT): Wind turbine fusedwind object

        Outputs:
           WindTurbine (WindTurbine)
        """
        self.name = name
        self.H = wt_desc.hub_height
        self.R = wt_desc.rotor_diameter / 2.0
        self.CT_idle = CT_idle

        #refCurvesArray=np.loadtxt(refCurvesFile,delimiter=', ',skiprows=5)
        #self.refCurvesArray = refCurvesArray
        # sanity check on the power curve and c_t_curve resamblance
        np.testing.assert_almost_equal(wt_desc.power_curve[:,0], wt_desc.c_t_curve[:,0])
        self.refCurvesArray = np.vstack([wt_desc.power_curve[:,0], wt_desc.power_curve[:,1], wt_desc.c_t_curve[:,1]]).T

        self.wt_init()

        #self.ref_u = wt_desc.power_curve[:,0]
        #self.ref_P = wt_desc.power_curve[:,1]
        #self.ref_CT = wt_desc.c_t_curve[:,1]

        #self.P_rated = np.max(self.ref_P)
        #self.u_cutin = wt_desc.cut_in_wind_speed
        #self.u_cutout = wt_desc.cut_out_wind_speed

def generate_GenericWindTurbinePowerCurveVT(WT):
    """Generate a GenericWindTurbinePowerCurveVT instance from a WindTurbine
    instance

    Parameters
    ----------
    WT: WindTurbine
        A WindTurbine instance

    Returns
    -------
    wt_desc: GenericWindTurbinePowerCurveVT
        A GenericWindTurbinePowerCurveVT instance containing the same information
        as in WT

    """
    wt_desc = GenericWindTurbinePowerCurveVT()
    wt_desc.hub_height = WT.H
    wt_desc.rotor_diameter = WT.R*2.0
    wt_desc.c_t_curve = np.vstack([WT.ref_u, WT.ref_CT]).T
    wt_desc.power_curve = np.vstack([WT.ref_u, WT.ref_P]).T
    wt_desc.cut_in_wind_speed = WT.u_cutin
    wt_desc.cut_out_wind_speed = WT.u_cutout
    return wt_desc

def generate_WTPC(wt, name, position):
    """Generate a WTPC instance from a WindTurbine
    instance

    Parameters
    ----------
    WT: WindTurbine
        A WindTurbine instance

    Returns
    -------
    wt_desc: WTPC
        A GenericWindTurbinePowerCurveVT instance containing the same information
        as in WT

    """
    wt_desc = WTPC()
    wt_desc.hub_height = wt.H
    wt_desc.rotor_diameter = wt.R*2.0
    wt_desc.c_t_curve = np.vstack([wt.ref_u, wt.ref_CT]).T
    wt_desc.power_curve = np.vstack([wt.ref_u, wt.ref_P]).T
    wt_desc.cut_in_wind_speed = wt.u_cutin
    wt_desc.cut_out_wind_speed = wt.u_cutout
    wt_desc.name = name
    wt_desc.position = position
    return wt_desc

class FWindFarm(WindFarm):
    def __init__(self, name, wt_layout):
        """Initializes a WindFarm object

        Inputs:
           name (str):  WindFarm name
           wt_layout (GenericWindFarmTurbineLayout): Wind Farm layout object

        Outputs:
           WindTurbine (WindFarm)
        """
        self.name = name
        self.pos = wt_layout.wt_positions.T
        self.nWT = self.pos.shape[1]
        assert_equal(self.pos.shape, [2, wt_layout.n_wt], err_msg='self.pos.shape should be [2, nWT]')


        # Vector from iWT to jWT: self.vectWTtoWT[:,i,j]
        self.vectWTtoWT=np.swapaxes([self.pos - np.repeat(
                np.atleast_2d(self.pos[:,i]).T, self.nWT, axis=1) for i in range(self.nWT)], 0, 1)


        # self.Xpos = coordArray[:,0]
        # self.Ypos = coordArray[:,1]
        # self.nWT = len(self.Ypos)
        # # TODO: Figure out a better way to add the turbine name
        self.WT = FWindTurbine('wt', wt_layout.wt_list[0])

def generate_GenericWindFarmTurbineLayout(WF):
    """Generate a GenericWindFarmTurbineLayout instance based on a WindFarm instance

    Parameters
    ----------
    WF: WindFarm
        A WindFarm instance

    Returns
    -------
    wt_layout: GenericWindFarmTurbineLayout
        A GenericWindFarmTurbineLayout instance
    """
    wt_layout = GenericWindFarmTurbineLayout([generate_WTPC(wt=WF.WT, position=WF.pos[:,i], name='wt%d'%(i+1)) for i in range(WF.nWT)])
    return wt_layout

# Rosetta stone to convert pure python inputs into FUSED-Wind wrapper inputs
rosetta = {
    'WS' : 'wind_speed',
    'z0' : 'roughness',
    'TI' : 'TI',
    'WD' : 'wind_direction',
    'NG' : 'NG',
    'sup' : 'sup',
    'pars' : 'pars'}

@implement_base(GenericWindFarm)
class FGCLarsen(Component):
    """Implementation of the Gunner C. Larsen stationary wake model in according
    to fusedwind.plant_flow interface.
    """

    # Interface from GenericWindFarm
    # Inputs:
    wind_speed = Float(iotype='in', low=0.0, high=100.0, units='m/s',
        desc='Inflow wind speed at hub height')
    wind_direction = Float(iotype='in', low=0.0, high=360.0, units='deg',
        desc='Inflow wind direction at hub height')
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='in',
        desc='wind turbine properties and layout')

    # Outputs:
    power = Float(iotype='out', units='kW',
        desc='Total wind farm power production')
    thrust = Float(iotype='out', units='N',
        desc='Total wind farm thrust')
    wt_power = Array([], iotype='out',
        desc='The power production of each wind turbine')
    wt_thrust = Array([], iotype='out',
        desc='The thrust of each wind turbine')

    # Specific to GCLarsen
    # Inputs:
    pars = List([0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0], iotype='in',
        desc='GCLarsen model parameters')
    TI = Float(0.05, iotype='in', low='0.0', high='1.0',
        desc='Ambient turbulence intensity')
    roughness = Float(0.0001, iotype='in', units='m', low='0.0', high='1.0',
        desc='The roughness height')
    NG = Int(4, low=4, high=6, iotype='in',
        desc='number of points in the Gauss integration')
    sup = Enum('lin', ['lin', 'quad'], iotype='in',
        desc='Wake supperposition model linear/quadratic')
    windfarm_name = Str('windfarm', iotype='in',
        desc='The wind farm name')
    windturbine_name = Str('windturbine', iotype='in',
        desc='The wind turbine name')

    # Outputs:
    wt_wind_speed = Array([], iotype='out', units='m/s',
        desc='The equivalent hub wind speed at each wind turbine')

    def execute(self):
        """ do your calculations here """
        # Create the WF object from the wt_layout
        WF = FWindFarm(self.windfarm_name, self.wt_layout)

        # Run the wind case
        self.wt_power, self.wt_wind_speed, wt_CT  = gcl.GCLarsen(
            WF = WF,
            WS = self.wind_speed,
            z0 = self.roughness,
            TI = self.TI,
            WD = self.wind_direction,
            NG = self.NG,
            sup = self.sup,
            pars = self.pars)

        self.wt_thrust = wt_CT * 0.5 * 1.225 * self.wt_wind_speed**2.0 * np.pi * WF.WT.R**2.

        self.power = self.wt_power.sum()
        self.thrust = self.wt_thrust.sum()
