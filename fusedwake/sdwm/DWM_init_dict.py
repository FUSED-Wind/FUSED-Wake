# -*- coding: utf-8 -*-
""" Initialization of storage dictionnaries from wind farm model and core DWM model
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""
import numpy as np

def init(WF):
    """Function that initializes dictionnaries used as part of the wake tree model (windfarm model) and the core DWM model

        Parameters
        ----------
        WF: WindFarm
            Instance of class WindFarm holding wind farm layout and coordinate

        Returns
        -------
        deficits: dict(nWT)
            holding a list of deficits contributions from upstream wakes
        turb:  dict(nWT)
            holding a list of turbulence intensities contributions from upstream
            wakes
        inlets_ffor: dict(nWT)
            holding a list of array containing the flow field in the fixed frame
            of reference from upstream wakes contributions
        inlets_ffor_deficits: dict(nWT)
            holding a list of array containing the flow field in the fixed frame
            of reference from upstream wakes contributions at the rotor position
        inlets_ffor_turb: dict(nWT)
            holding a list of array containing the turbulence field in the fixed
            frame of reference from upstream wakes contributions at the rotor
            position
        out: dict(nWT)
            holding the main sDWM model outputs (i.e mean power from BEM, mean
            power estimated from power curve, mean rotor averaged wind speed,
            mean rotor average turbulence intensity, mean thrust coefficient
            from BEM and from power curve)
        DWM: dict(nWT)
            list containing full outputs of the sDWM (including flow field in
            ffor and mfor) See description of DWM_outputs for more details
        ID_waked: dict(nWT)
            holding list of upstream turbine index for each turbine in the wind
            farm
        ID_wake_adj: dict(nWT)
            holding list of downstream turbine index for each turbine in the wind
            farm
        Farm_p_out: float
            wind farm power produced
        WT_p_out: np.array(float)
            power produced by each turbine in the wind farm
        Vel_out: np.array(float)
            mean rotor averaged aggregated velocity for each turbine in the wind
            farm
        Pos_out: np.array(float)
            outputs coordinates X,Y of wind turbines
    """
    ## Init dictionnaries involved in full outputs
    deficits = dict()
    for iC in range(WF.nWT):
        deficits[str(iC)]=[]  # Deficits tree dictionnary init
    turb = dict()
    for iC in range(WF.nWT):
        turb[str(iC)]=[]    ## Turbulence intensities tree dictionnary init
    inlets_ffor = dict() ## inlets in fixed frame of reference to the wake model
    for iC in range(WF.nWT):
        inlets_ffor[str(iC)]=[]
    inlets_ffor_deficits = dict()  ## deficits in fixed frame of reference to the wake model
    for iC in range(WF.nWT):
        inlets_ffor_deficits[str(iC)]=[]
    inlets_ffor_turb = dict()
    for iC in range(WF.nWT):
        inlets_ffor_turb[str(iC)]=[]    ## inlets dictionnary for turbulence build up
    out = dict()    ## outputs dictionnary to the wake model
    DWM = dict()## flow field dictionnary outputs to wake model
    ID_waked=dict()
    for iC in range(WF.nWT):
        ID_waked[str(iC)]=[]    # ID_Wake ordered by wind direction
    ID_wake_adj=dict()
    for iC in range(WF.nWT):
        ID_wake_adj[str(iC)]=[]


    # Power output list
    Farm_p_out=0.
    WT_p_out=np.zeros((WF.nWT))
    Vel_out=np.zeros((WF.nWT))
    Pos_out=np.zeros((WF.nWT,2))
    return deficits, turb, inlets_ffor, inlets_ffor_deficits,inlets_ffor_turb,out, DWM, ID_waked, ID_wake_adj, Farm_p_out, WT_p_out, Vel_out, Pos_out
