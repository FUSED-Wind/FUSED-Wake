# -*- coding: utf-8 -*-
""" Misc functions from variable sources
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import pandas as pd


def LoadOutputs(folder,vWD,WF,WS,TI):
    powers=np.zeros((len(range(WF.nWT)),len(vWD)))
    ref_powers=[]
# for iD in range(0,1):
    for iD in np.arange(0,len(vWD),1):
        WD=float(vWD[iD])
        filename= folder +'/dwm_WS_'+ str(WS) +'_dir_' + str(WD)+'_TI_' + str(TI) +'.npy'
        # print filename
        try:
            out=np.load(filename).item()
            for iK in range(WF.nWT):
                # powers[iK,iD]=out[str(iK)][0] # BEM
                powers[iK,iD]=out[str(iK)][4] # power curve
            ref_powers.append(max(powers[:,iD]))
        except:
            for iK in range(WF.nWT):
                powers[iK,iD]=np.nan
            ref_powers.append(np.nan)

    return powers, ref_powers


# Pierre Elouan Rethore functions
def my_rolling_deg(df, x='wd', y='eff', dwd=2.5):
    inte = interp1d(df[x], df[y])
    inte2 = lambda x_: inte(ism360(x_, df[x].max()))
    def filter_func(d):
        return {y:quad(inte2, d-dwd, d+dwd)[0]/(2.*dwd),x:d}
    return pd.DataFrame(map(filter_func, df[x]))


# Pierre Elouan Rethore functions
def ism360(v, endp):
    """Make sure the direction is in [0.0, 360]"""
    if np.isnan(v):
        return v
    if v>=0.0:
        if v>endp and v<360.0:
            return endp
        elif v>=360.0:
            return v - 360.0
        else:
            return v
    else:
        return ism360(v + 360.0, endp)




def smooth(x,window_len=11,window='hanning'):
        """
        Perform smoothing (moving average) from input vector
        """
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


def to_bool(value):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))
