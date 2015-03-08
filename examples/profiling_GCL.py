from numpy import linspace
from fusedwind.plant_flow.asym import AEPMultipleWindRoses
from gclarsen.fusedwasp import PlantFromWWH
from gclarsen.fused import FGCLarsen
import cProfile

__author__ = 'pire'

def prof():
    ### Single wind rose type
    hrAEP = AEPMultipleWindRoses()
    hrAEP.add('wf', FGCLarsen())
    hrAEP.configure()
    #hrAEP.connect('wt_layout', 'wf.wt_layout')
    wt_layout = PlantFromWWH(filename='hornsrev1_turbine_nodescription.wwh').wt_layout
    hrAEP.wf.wt_layout = wt_layout  ## <- to speed things up
    hrAEP.wt_layout = wt_layout
    hrAEP.wind_directions = linspace(0.0, 360.0, 36)[:-1]
    hrAEP.wind_speeds =linspace(4.0, 25.0, 22)
    hrAEP.run()
    print hrAEP.net_aep
    print hrAEP.wt_aep
cProfile.run('print prof(); print', 'restats')

import pstats
p = pstats.Stats('restats')
p.sort_stats('cumulative').print_stats(100)
print 'done'