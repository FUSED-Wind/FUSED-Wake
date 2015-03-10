from numpy import hstack
from openmdao.main.component import Component
from openmdao.main.datatypes.array import Array
from openmdao.main.datatypes.str import Str
from openmdao.main.datatypes.vtree import VarTree
from fusedwind.interface import implement_base
from fusedwind.plant_flow.asym import BaseAEPModel, AEPWindRose
from fusedwind.plant_flow.vt import WTPC, GenericWindRoseVT, GenericWindTurbinePowerCurveVT, \
    GenericWindFarmTurbineLayout

__author__ = 'pire'


from py4we.wasp import WWF, WTG, WWH


class WTDescFromWTG(Component):

    """Create a wt_desc from a .wtg WAsP file"""
    # Inputs
    filename = Str(iotype='in',
        desc='The .wtg file name')

    # Outputs
    wt_desc = VarTree(GenericWindTurbinePowerCurveVT(), iotype='out',
        desc='The wind turbine power curve')


    def __init__(self, filename=None):
        super(WTDescFromWTG, self).__init__()

        if filename is not None:
            self.filename = filename
            self.execute()


    def execute(self):
        # Reading the .wtg file
        wtg = WTG(self.filename)

        # Filling up the VT
        self.wt_desc.power_curve = wtg.data[:, :2]
        self.wt_desc.c_t_curve = hstack([wtg.data[:, 0], wtg.data[:, 2]])
        self.wt_desc.cut_in_wind_speed = wtg.data[0, 0]
        self.wt_desc.cut_out_wind_speed = wtg.data[-1, 0]
        self.wt_desc.air_density = wtg.density
        self.wt_desc.rotor_diameter = wtg.rotor_diameter
        self.wt_desc.hub_height = wtg.hub_height
        self.wt_desc.test_consistency()


class PlantFromWWF(Component):

    """Create a Plant information from a .wwf WAsP file"""
    # Inputs
    filename = Str(iotype='in',
        desc='The .wwf file name')
    wt_desc = VarTree(GenericWindTurbinePowerCurveVT(), iotype='in',
        desc='The wind turbine power curve')

    # Outputs
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='out',
        desc='wind turbine properties and layout')
    wind_rose_array = Array(iotype='out', units='m/s',
        desc='Windrose array [wind_directions, frequency, weibull_A, weibull_k]')

    def execute(self):
        # Reading the .wwf file
        wwf = WWF(self.filename)
        self.wt_layout.wt_list.append(self.wt_desc)
        # self.wt_layout.wt_wind_roses.frequency_array =
        self.wt_layout.configure_single()

        for wt, wr in self.wwf.windroses.iteritems():
            self.wt_layout.wt_positions[i, :] = self.wwf.data[wt][:2]
            self.wt_layout.wt_wind_roses.append(wr)
            i += 1

        self.wt_layout.configure_single()
        # For practical reason we also output the first wind rose array outside
        # the wt_layout
        self.wind_rose_array = self.wt_layout.wt_wind_roses[0]


class PlantFromWWH(Component):

    """Create a Plant information from a .wwh WAsP file"""
    # Inputs
    filename = Str(iotype='in',
        desc='The .wwh file name')

    # Outputs
    wt_layout = VarTree(GenericWindFarmTurbineLayout(), iotype='out',
        desc='wind turbine properties and layout')
    # wind_rose_array = Array(iotype='out', units='m/s',
    #     desc='Windrose array [wind_directions, frequency, weibull_A, weibull_k]')

    def __init__(self, filename=None):
        """
        The constructor that can take the wwh filename as optional parameter.
        :param filename: wwh file name [optional]
        :return: self
        """
        super(PlantFromWWH, self).__init__()
        if filename is not None:
            self.filename = filename
            self.execute()

    def execute(self):
        # Reading the .wwf file
        wwh = WWH(self.filename)

        for wt, data in wwh.wind_turbines.iteritems():
            turbine = wwh.turbine_descriptions[data['type']]
            self.wt_layout.add_wt(WTPC(
                name='wt_'+wt,
                position=data['position'][:2],
                wind_rose=GenericWindRoseVT(weibull_array=data['wind_rose']),
                power_curve = turbine['data'][:, :2],
                c_t_curve = turbine['data'][:, [0, 2]],
                cut_in_wind_speed = turbine['data'][0, 0],
                cut_out_wind_speed = turbine['data'][-1, 0],
                air_density = turbine['density'],
                rotor_diameter = turbine['rotor_diameter'],
                hub_height = data['position'][2]
            ))


# @implement_base(BaseAEPModel, AEPWindRose)
# class WWHAEP(FUSEDAssembly):
#
#     """Class that loads a wind farm position and wind roses from a .wwh WAsP file
#     and perform an AEP calculation.
#     """
#
#     wind_rose_type = Enum('single', ('single', 'multiple'),
#         desc='Are we using only one wind rose for the whole wind farm, or a different wind rose for each turbine?')
#
#     # Slots
#     wwh = Slot(PlantFromWWH)
#
#     # Interface Slots (using the @implement_base)
#     wf = InterfaceSlot(GenericWindFarm,
#         desc='A wind farm assembly or component')
#     postprocess_wind_rose = InterfaceSlot(GenericPostProcessWindRose,
#         desc='The component taking care of postprocessing the wind rose')
#     case_gen = InterfaceSlot(GenericWindRoseCaseGenerator,
#         desc='Generate the cases from the inputs')
#
#     # Inputs
#     filename = Str(iotype='in',
#         desc='The .wwh file name')
#     wind_speeds = List([], iotype='in', units='m/s',
#         desc='The different wind speeds to run [nWS]')
#     wind_directions = List([], iotype='in', units='deg',
#         desc='The different wind directions to run [nWD]')
#
#     # Outputs
#     array_aep = Array([], iotype='out', units='kW*h',
#         desc='The energy production per sector [nWD, nWS]')
#     gross_aep = Float(iotype='out', units='kW*h',
#         desc='Gross Annual Energy Production before availability and loss impacts')
#     net_aep = Float(iotype='out', units='kW*h',
#         desc='Net Annual Energy Production after availability and loss impacts')
#     capacity_factor = Float(0.0, iotype='out',
#         desc='Capacity factor for wind plant')
#
#     def configure(self):
#         # Adding
#         self.add('wwh', PlantFromWWH())
#         self.add('')
#
#         # Configure the wind rose postprocessor: Are we using only one wind rose
#         # for the whole wind farm, or a different wind rose for each turbine?
#         if self.wind_rose_type == 'single':
#             self.add('case_gen', SingleWindRoseCaseGenerator())
#             self.add('postprocess_wind_rose', PostProcessSingleWindRose())
#         if self.wind_rose_type == 'multiple':
#             self.add('postprocess_wind_rose', PostProcessMultipleWindRoses())
#             self.add('case_gen', MultipleWindRosesCaseGenerator())
#
#         # Base class configure
#         configure_AEPWindRose(self)
#
#         # Wiring
#         self.connect('wwh.wt_layout', 'wf.wt_layout')
#         self.connect('filename', 'wwh.filename')
#
#         if self.wind_rose_type == 'single':
#             self.connect('wwh.wind_rose', 'case_gen.wind_rose')
#         if self.wind_rose_type == 'multiple':
#             self.connect('wwh.wt_layout', 'case_gen.wt_layout')
#
#         # Re-organizing the workflow
#         self.driver.workflow.add(['wwh', 'case_generator', 'wind_rose_driver',
#                                   'postprocess_wind_rose'])

#
#
# @implement_base(BaseAEPModel, AEPWindRose)
# class WWHAEP(FUSEDAssembly):
#
#     """Class that loads a wind farm position and wind roses from a .wwh WAsP file
#     and perform an AEP calculation.
#     """
#
#     wind_rose_type = Enum('single', ('single', 'multiple'),
#         desc='Are we using only one wind rose for the whole wind farm, or a different wind rose for each turbine?')
#
#     # Slots
#     wwh = Slot(PlantFromWWH)
#
#     # Interface Slots (using the @implement_base)
#     wf = InterfaceSlot(GenericWindFarm,
#         desc='A wind farm assembly or component')
#     postprocess_wind_rose = InterfaceSlot(GenericPostProcessWindRose,
#         desc='The component taking care of postprocessing the wind rose')
#     case_gen = InterfaceSlot(GenericWindRoseCaseGenerator,
#         desc='Generate the cases from the inputs')
#
#     # Inputs
#     filename = Str(iotype='in',
#         desc='The .wwh file name')
#     wind_speeds = List([], iotype='in', units='m/s',
#         desc='The different wind speeds to run [nWS]')
#     wind_directions = List([], iotype='in', units='deg',
#         desc='The different wind directions to run [nWD]')
#
#     # Outputs
#     array_aep = Array([], iotype='out', units='kW*h',
#         desc='The energy production per sector [nWD, nWS]')
#     gross_aep = Float(iotype='out', units='kW*h',
#         desc='Gross Annual Energy Production before availability and loss impacts')
#     net_aep = Float(iotype='out', units='kW*h',
#         desc='Net Annual Energy Production after availability and loss impacts')
#     capacity_factor = Float(0.0, iotype='out',
#         desc='Capacity factor for wind plant')
#
#     def configure(self):
#         # Adding
#         self.add('wwh', PlantFromWWH())
#
#         # Configure the wind rose postprocessor: Are we using only one wind rose
#         # for the whole wind farm, or a different wind rose for each turbine?
#         if self.wind_rose_type == 'single':
#             self.add('case_gen', SingleWindRoseCaseGenerator())
#             self.add('postprocess_wind_rose', PostProcessSingleWindRose())
#         if self.wind_rose_type == 'multiple':
#             self.add('postprocess_wind_rose', PostProcessMultipleWindRoses())
#             self.add('case_gen', MultipleWindRosesCaseGenerator())
#
#         # Base class configure
#         configure_AEPWindRose(self)
#
#         # Wiring
#         self.connect('wwh.wt_layout', 'wf.wt_layout')
#         self.connect('filename', 'wwh.filename')
#
#         if self.wind_rose_type == 'single':
#             self.connect('wwh.wind_rose', 'case_gen.wind_rose')
#         if self.wind_rose_type == 'multiple':
#             self.connect('wwh.wt_layout', 'case_gen.wt_layout')
#
#         # Re-organizing the workflow
#         self.driver.workflow.add(['wwh', 'case_generator', 'wind_rose_driver',
#                                   'postprocess_wind_rose'])
