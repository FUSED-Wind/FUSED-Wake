import numpy as np
import sys
sys.path.append("../src/")
import numpy as np
import unittest

from fusedwake.sdwm.DWM_flowfield_farm import run_sdwm


class Test_Lillgrund_SDWM(unittest.TestCase):

    def test_120(self):
        np.random.seed(1) # make sure that the turbulence is the same
        inputs={'WD':120.0,
            'WS':9.0,
            'TI':0.06,
            # WTcoord='../data/Lill_A04B04C04.dat'; '../data/Lill_full.dat'
            # 'WTcoord':'../data/Lill_A04.dat',
            'WTcoord':'../data/Lill_full.dat',
            'WTG':'NREL5MW',
            'HH':90.0,
            'R':63.0,
            'stab':'N',
            'accum':'dominant'}


        result_120 = np.array([2696.42,  2696.42,  2696.42,  2696.42,   597.58,  2696.42,   414.54,  2696.42,
                               423.7,    2696.42,  443.03,   470.26,   419.33,   379.35,   448.09,   410.47,
                               417.32,   382.25,   2696.42,  483.34,   415.44,   481.11,   391.84,   445.05,
                               406.9,    484.04,   417.47,   702.95,   497.09,   483.47,   538.81,   449.31,
                               503.19,   378.21,   577.27,   504.95,   519.57,   499.66,   569.63,   1628.06,
                               591.5,    542.,     588.37,   593.08,   460.87,   644.2,    605.99,   613.35,])

        Farm_p_out,WT_p_out,Vel_out,Pos_out,WF,WT,aero, meta, mfor, ffor, DWM, deficits,inlets_ffor,inlets_ffor_deficits, inlets_ffor_turb,turb, out,ID_waked = run_sdwm(**inputs)
        np.testing.assert_array_almost_equal(np.array(WT_p_out), result_120)


if __name__ == "__main__":
    unittest.main()
