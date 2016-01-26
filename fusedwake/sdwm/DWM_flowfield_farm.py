# -*- coding: utf-8 -*-
""" Main DWM core program: flowfield calculation

@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
from cDWM import Meta, Aero, Meand, FFoR, MFoR, Outputs
from cBEM import InitBEM
from DWM_calc_mixL import DWM_calc_mixL
from DWM_main_BEM import getInduction
from math import pi, sqrt
from scipy import io, interpolate
from DWM_misc import smooth

def DWM_main_field_model(ID_waked,deficits,inlets_ffor,inlets_ffor_deficits,inlets_ffor_turb,turb,DWM,out,**par):
    """Main flow field calculation function, handling all the calls to each sub functions. This function is called for
       each turbine in the wind farm from the most upstream to the most downstream one. The flow field calculations are
       carried out at the downstream distance of interest, i.e., where downstream rotors are partially or fully in the
       wake of the upstream turbine.


    Parameters
    ----------
    ID_waked: dict(nWT)
        holding list of upstream turbine index for each turbine in the wind farm
    deficits: dict(nWT)
        holding a list of the mean rotor averaged deficit from upstream wakes
    turb: dict(nWT)
        holding a list of mean turbulence intensities contributions from upstream wakes
    inlets_ffor: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position
    inlets_ffor_turb: dict(nWT)
        holding a list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position
    out: dict(nWT)
        holding the main sDWM model outputs i.e mean power from BEM, mean power estimated from powercurve, mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and from power curve
    par: dict(nWT)
        holding the DWM cylindrical and cartesian grid coordinates as well as turbine and ambient conditions


    Returns
    -------
    meta: (instance of class)
        Instance of class Meta holding DWM core variables
    aero: (instance of class)
        Instance of class Aero holding BEM-aero core variables
    mfor: (instance of class)
        Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model in local Mixl coordinates
    ffor: (instance of class)
        Instance of class Ffor holding the fixed frame of reference velocity field in global WF coordinates
    DWM: dict(nWT)
        list containing full outputs of the sDWM (including flow field in ffor and mfor) See description of DWM_outputs for more details
    deficits: dict(nWT)
        update list of deficits contributions from upstream wakes
    inlets_ffor: dict(nWT)
        updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT)
        updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position
    inlets_ffor_turb: dict(nWT)
        updated list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position
    turb: dict(nWT)
        updated list of mean turbulence intensities contributions from upstream wakes
    out: dict(nWT)
        returns the mean power from BEM, mean power estimated from power curve,mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and from power curve
    ID_waked: dict(nWT)
        holding list of upstream turbine index for each turbine in the wind farm
    """

    ## Create instances of class
    meta=Meta()
    meta.parse(**par)
    meand=Meand()
    ffor =  FFoR()
    aero = Aero(meta.WTG)
    ###### Set up MFoR and FFoR streamwise domain properties   #########################################################
    meta                 = DWM_make_grid(meta)
    ###### Load wake meandering properties from meta model: f(stab,hub height,z,TI) ####################################
    meand                = DWM_meta_meand(meand,meta)
    ###  Run BEM model and create velocity deficit calculations inputs #################################################
    aero,mfor,out,BEM = DWM_aero(meta,ffor,aero,deficits,turb,inlets_ffor,inlets_ffor_deficits,out)
    ############### Perform wake velocity calculations in MFoR #########################################################
    mfor                 = DWM_calc_mixL(meta,aero,mfor)
    ############## Reconstruct global flow field by applying wake meandering ###########################################
    ffor,meta            = DWM_MFOR_to_FFOR(mfor,meta,meand,ffor)
    ############### Compute deficit at downstream rotor ################################################################
    deficits, ID_waked,inlets_ffor,inlets_ffor_deficits  \
                         = DWM_get_deficit(ffor,meta,deficits,ID_waked,inlets_ffor,inlets_ffor_deficits)
    ############### Compute turbulence level at downstream rotor #######################################################
    turb,inlets_ffor_turb                 = DWM_get_turb(ffor,meta,turb,inlets_ffor_turb)
    ############## Write relevant results in DWM variables #############################################################
    if meta.full_output is True:
        DWM                  = DWM_outputs(DWM,ffor,mfor,meta,aero,BEM)
    return( aero, meta, mfor, ffor, DWM, deficits, inlets_ffor, inlets_ffor_deficits, inlets_ffor_turb, turb, out, ID_waked)

def DWM_aero(meta,ffor,aero,deficits,turb,inlets_ffor,inlets_ffor_deficits,out):
    """ Aerodynamique module of the DWM. This module contains the wake summation module (deficit and turbulence accumulation) as well as a steady state blade element momentum

    Parameters
    ----------
    meta: (instance of class)
        Instance of class Meta holding DWM core variables
    aero: (instance of class)
        Instance of class Aero holding BEM-aero core variables

    deficits: dict(nWT)
        holding a list of deficits contributions from upstream wakes
    turb:  dict(nWT)
        holding a list of mean turbulence intensities contributions from upstream wakes
    inlets_ffor: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of reference from
        upstream wakes contributions
    inlets_ffor_deficits: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
    out: dict(nWT)
        holding the main sDWM model outputs i.e mean power from BEM, mean power estimated from powercurve,
        mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and
        from power curve

    Returns
    -------
    aero: (instance of class)
        updated Instance of class Aero holding BEM-aero core variables
    mfor: (instance of class)
        updated Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
    out dict(nWT):
        dict including mean power from PC and BEM, mean thrust coefficient from PC and BEM
    BEM: (instance of class)
        holds the key results from the BEM calculation

    """

    mfor   =  MFoR(meta.WTG)
    ## Compute the average wake deficit accumulation
    if meta.wake_ind_setting==1:
        rWS = deficits.get(str(meta.wtg_ind[0]))
        if not rWS: # means that the current turbine is in free stream
            rWS=np.array([1.0])
        if meta.accu == 'linear':
            meta.mean_WS_DWM=meta.WS*(1.-(np.sum([1. - xx for xx in rWS])))
        elif meta.accu == 'quadratic':
            meta.mean_WS_DWM=meta.WS*(1.-np.sqrt(np.sum([(1.-xx)**2 for xx in rWS])))
        elif meta.accu == 'dominant':
            meta.mean_WS_DWM= meta.WS*(1.-(np.max([1. - xx for xx in rWS])))
        else:
            print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
    else: #  'bypassed':
        meta.mean_WS_DWM   = meta.WS # free stream everywhere similar to current implementation of HAWC2

    # Set buildup of turbulence
    if meta.Tbuildup_setting==1:
        ti = turb.get(str(meta.wtg_ind[0]))
        if not ti:
            meta.mean_TI_DWM  = meta.TI
        else:
            meta.mean_TI_DWM  = np.max(ti)
    else:
        meta.mean_TI_DWM  = meta.TI

    # Run BEM at accumulated deficit
    aero,BEM   =  DWM_rotor_aero(meta,aero)

    # domain induction
    a_domain     = np.interp(meta.vr_m,np.hstack(([aero.r_w, aero.r_w[-1]+0.01, aero.r_w[-1]+0.02])), np.hstack((( 1.0 - aero.U_w), [0., 0.])))

    ## Compute the accumulated flow field for accurate inlet definition of the MFoR wake calculation
    # center all disks before wake accumulation
    if not inlets_ffor.get(str(meta.wtg_ind[0])):# upstream rotor
        # mfor.U_init =
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    elif meta.accu_inlet is False: # HAWC2-DWM behavior
        # mfor.U_init =   None
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    else:   # if turbine not in the freestream, we need to compute the proper accumulated inlet to the turbine
        ranger=np.linspace(-1.,1.,meta.dR*2.)  # np.linspace(-2.,2.,meta.dR*4.)
        inlets_ffor_deficits_np_3D=np.ones((len(ranger) ,len(ranger)   , len(inlets_ffor[str(meta.wtg_ind[0])])))
        grid_x, grid_y = np.mgrid[-1.:1.:meta.dR*2j, -1.:1.:meta.dR*2j]
        for ii in range(len(inlets_ffor[str(meta.wtg_ind[0])])):
            offsets=2.*(min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])+abs(max(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])-min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]))/2.)
            # need to interp on a new array of equal size
            values=inlets_ffor_deficits[str(meta.wtg_ind[0])][ii]
            X, Y = np.meshgrid(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]-offsets,inlets_ffor[str(meta.wtg_ind[0])][ii][0][1])
            points=np.vstack((np.ravel(X),np.ravel(Y)))
            wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='linear')
            inlets_ffor_deficits_np_3D[:,:,ii]=wake_i
        # flow field superposition
        if meta.accu == 'linear':
            U_init=1.-(np.sum(1.-inlets_ffor_deficits_np_3D,axis=2))
        elif meta.accu == 'quadratic':
            U_init=1.-np.sqrt(np.sum((1.-inlets_ffor_deficits_np_3D)**2,axis=2))
        elif meta.accu == 'dominant':
            U_init=np.amin(inlets_ffor_deficits_np_3D, axis=2)
        else:
            print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
        # Transform to axisymmetric profile inlet, required by Ainslie model (N-S TL model)
        r_dist_2= np.sqrt(grid_x**2 + grid_y**2 )  #Find distance to centre of wake plane
        ffor.WS_axial_sym      = np.ones((len(np.arange(0,meta.dR+1.,1))))
        ffor.WS_axial_sym[0]=np.nanmean(U_init[r_dist_2 < (1.05*np.amin(r_dist_2))])
        for i_r_pos in np.arange(1,meta.dR+1.,1):
            a=r_dist_2 > ((i_r_pos+1-1.5)*(1.0/meta.dR))# rotor neg boundaries
            bb=r_dist_2 < ((i_r_pos+1-0.5)*(1.0/meta.dR)) #rotor pos boundaries
            c=np.logical_and(a,bb)
            bin_filter = c
            tmp_ffor_flow_field_ws_mean          = U_init[bin_filter]
            ffor.WS_axial_sym[i_r_pos]           = np.nanmean(tmp_ffor_flow_field_ws_mean)
        ffor.r_sym=np.arange(0,meta.dR+1.,1)/meta.dR
        # Update the DWM inlet
        if ffor.r_sym[-1] >= meta.vr_mixl[-1]:
            # mfor.U_init = (1.0-a_domain) * np.interp(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
            mfor.U_init = (1.0-a_domain) * interpolate.InterpolatedUnivariateSpline(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
        else:
            mfor.U_init = (1.0-a_domain) * np.hstack((ffor.WS_axial_sym.ravel(), np.ones((((meta.dR * meta.lr_mixl)-ffor.WS_axial_sym.size),1)).ravel()))
        # Finishing
        mfor.U_init[mfor.U_init < 0.0]=0.0 # prevent from negative velocities on linear summation
        mfor.U_init=smooth( mfor.U_init,window_len=5)

    # Power curve based
    try:
        aero.pow_cur=meta.WTG_spec.get_P(meta.mean_WS_DWM)
        aero.ct_cur=meta.WTG_spec.get_CT(meta.mean_WS_DWM)
    except:
        aero.pow_cur=0.
        aero.ct_cur=0.
    # write outlets
    out[str(meta.wtg_ind[0])]=[]
    out[str(meta.wtg_ind[0])].append(float(format(aero.Power/1000., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_WS_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_TI_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.CT/1., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.pow_cur, '.2f'))) # based on poewr curve
    out[str(meta.wtg_ind[0])].append(float(format(aero.ct_cur, '.2f'))) # based on poewr curve
    return aero, mfor, out, BEM

def DWM_rotor_aero(meta,aero):
    """
    Function that calls the BEM calculation and calculate the corrected initial wake radius
    (boundary condition to Ainslie model ref [2])

    Parameters
    ----------
    aero:    class holding the aero parameters to the Ainslie model
    meta:    class holding grid and ambient parameters

    Returns
    -------
    aero (updated)
    CP [1],        : turbine power coefficient
    CPloc[r,1],    : local CP for power estimation
    Power [1],     : power from BEM [W]
    CT [1],        : thrust coefficient [-]
    RPM [1],       : rotational speed [RPM]
    PITCH [1]      : pitch angle [deg]
    U_w [1,r],     : local axial velocity in the wake from momemtum theory
    a [1,r],       : local induction factor
    dA [1],        : local annular area
    f_w [1],       : calibration factor on wake deficit expansion
    mean_a [1],    : rotor averaged axial induction
    r_t [1,r],     : non dimensional radial position where a is evaluated
    r_w [1,r]      : wake radius in near wake regime
    """

    if (meta.mean_WS_DWM >= meta.WTG_spec.u_cutin) or (meta.mean_WS_DWM <= meta.WTG_spec.u_cutout) is True:
        BEM = getInduction(30, meta.WTG, 'hawc', meta.mean_WS_DWM, meta) # 30 point per radius (hard coded value can be changed to input parsed param)
        aero.a = np.array(BEM.a)
        aero.r_t = np.array(BEM.r)/BEM.R
        aero.CP = np.array(BEM.CP)
        aero.Power = np.array(BEM.Power)
        aero.CT    = np.array(BEM.CT)
        aero.RPM    = BEM.RPM
        aero.PITCH    = BEM.PITCH

    else: # idle turbine simplified approach (not completely consistent)
        BEM = InitBEM(30)
        aero.r_t = np.arange(0.,1.+1./(meta.dR),1./(meta.dR))
        aero.Power  = 0.0 #kW
        aero.CT     = meta.WTG_spec.CT_idle
        aero.ia = meta.WTG_spec.get_a(aero.CT)
        aero.a=  aero.ia*np.ones(len(aero.r_t))

    # Boundary conditions for the r_w
    aero.dA = np.concatenate(([0], (pi*aero.r_t [1:]**2 - pi*aero.r_t [0:-1]**2)), axis=0)
    aero.mean_a= np.sum(aero.a*aero.dA)/pi
    # Uniform expansion
    aero.f_w        = sqrt( (1.0-aero.mean_a) / (1.0- ((1.0+meta.fR) * aero.mean_a))  )
    aero.r_w = np.dot(aero.r_t, aero.f_w)
    # Boundary conditions for U_w
    aero.U_w     = 1.0-(aero.a * (1.0 + meta.fU))

    return aero, BEM

def DWM_get_deficit(ffor,meta,deficits,ID_waked,inlets_ffor,inlets_ffor_deficits):
    """
    Function that performs the velocity integration for each downstream rotor of a given upstream turbine.

    Parameters
    ----------
    ffor:(instance of class)
        Instance of class Ffor holding the fixed frame of reference velocity field in global WF coordinates
    meta: (instance of class)
        Instance of class Meta holding DWM core variables
    deficits: dict(nWT)
        holding a list of deficits contributions from upstream wakes
    ID_waked: dict(nWT)
        holding list of upstream turbine index for each turbine in the wind farm
    inlets_ffor: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT)
        holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position

    Returns
    -------
    deficits: dict(nWT)
        updated list of deficits contributions from upstream wakes
    ID_waked: dict(nWT)
        updated list of upstream turbine index for each turbine in the wind farm
    inlets_ffor: dict(nWT)
        updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT)
        updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position
    """

    for i_z in np.arange(0,meta.nz,1):
        # on global frame mesh
        X,Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz=np.sqrt((X + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y)**2 )>=0.5
        wakedefmask = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=0.0).filled()
        wakedefmasknancoarse = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=np.nan).filled()
        disk = np.ma.array(np.zeros(wakedefmask.shape), mask=~index_trapz, fill_value=1.0).filled()
        disk_area=np.trapz(np.trapz(disk,dx=1./meta.dy),dx=1./meta.dx)
        trapz2=np.trapz(np.trapz(wakedefmask,dx=1./meta.dy),dx=1./meta.dx)

        deficits[str(meta.wtg_ind[i_z])].append(trapz2/disk_area)
        ID_waked[str(meta.wtg_ind[i_z])].append(meta.wtg_ind[0])
        inlets_ffor_deficits[str(meta.wtg_ind[i_z])].append(wakedefmasknancoarse)
        inlets_ffor[str(meta.wtg_ind[i_z])].append([np.vstack(((meta.x_vec-meta.hub_x[i_z]),(meta.y_vec-meta.hub_y),ffor.WS_axial_ffor[:,:,i_z]))])

    return deficits, ID_waked,inlets_ffor,inlets_ffor_deficits

def DWM_get_turb(ffor,meta,turb,inlets_ffor_turb,):
    """
    Function that calculate the rotor averaged turbulence intensity

    Parameters
    ----------
    ffor: (instance of class) Ffor holding the fixed frame of reference velocity field in global WF coordinates
    meta: (instance of class) Instance of class Meta holding DWM core variables
    turb:  dict(nWT) holding a list of turbulence intensities contributions from upstream wakes
    inlets_ffor_turb: dict(nWT) holding a list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position

    Returns
    -------
    turb:  dict(nWT)
        updated list of turbulence intensities contributions from upstream wakes
    inlets_ffor: dict(nWT)
        updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_turb: dict(nWT)
        updated list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position
    """

    for i_z in np.arange(0,meta.nz,1):
        X,Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz=np.sqrt((X + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y)**2 )>=0.5
        turbmask = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=0.0).filled()
        turbmasknan = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=np.nan).filled()
        disk = np.ma.array(np.zeros(turbmask.shape), mask=~index_trapz, fill_value=1.0).filled()
        disk_area=np.trapz( np.trapz(disk,dx=1./meta.dy),dx=1./meta.dx )
        trapz2=np.trapz( np.trapz(turbmask,dx=1./meta.dy),dx=1./meta.dx )
        turb[str(meta.wtg_ind[i_z])].append(trapz2/disk_area)
        inlets_ffor_turb[str(meta.wtg_ind[i_z])].append(turbmasknan)

    return turb,inlets_ffor_turb

def DWM_MFOR_to_FFOR(mfor,meta,meand,ffor):
    """
    Function that calculate the velocity in the fixed (global) frame of reference from the Mfor

    Parameters
    ----------
    meand: (instance of class)
        holding the meandering parameters
    meta: (instance of class)
        holding grid and ambient parameters

    Returns
    -------
    ffor: (instance of class)    updated class holding global velocity field
    x_vec(1,nx) : coordinates in X direction
    y_vec(1,ny) : coordinates in Y direction
    z_vec(1,nz) : coordinates in Z direction (streamwise)
    x_mat(nx,ny): coordinates matrix in meshgrid format for the X component
    y_mat(nx,ny): coordinates matrix in meshgrid format for the Ycomponent
    z_mat(nx,nz): coordinates matrix in meshgrid format for the Z component
    TI_meand_axial_ffor (nx,ny,nz): turbulence due to wake meandering in global coordinate system
    WS_axial_ffor (nx,ny,nz): velocity deficit in global coordinate system
    TI_axial_ffor (nx,ny,nz): apparent turbulence in global coordinate system see Madsen et al [2]
    meta: (instance of class)
    """

    ##############################################################################################################
    # recalculate into Cartesian grid
    # initiate/reset Cartesian flow field
    ffor.ffor_flow_field_TI_tmp_tmp =  meta.TI * np.ones((meta.nx,meta.ny))  #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.TI_axial_ffor_tmp     =  np.zeros((meta.nx,meta.ny,meta.nz)) # %X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.WS_axial_ffor_tmp     =  np.zeros((meta.nx,meta.ny,meta.nz)) # %X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.ffor_flow_field_ws_tmp2    =  np.zeros((meta.nx,meta.ny,meta.nz)) # %X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.TI_meand_axial_ffor    =  np.zeros((meta.nx,meta.ny,meta.nz))
    ffor.WS_axial_ffor =  np.zeros((meta.nx,meta.ny,meta.nz))
    ffor.TI_axial_ffor =  np.zeros((meta.nx,meta.ny,meta.nz))

    # CREATES THE GLOBAL FLOW FIELDS IN CARTESIAN GRID
    meta.z_mat = np.tile(meta.vz,(meta.nx,1))

    for i_z in np.arange(0,meta.nz,1):
        # EXTRACT TI_DWM AND WS_DWM IN MFoR
        try:
            DWM_WS_DATA = mfor.U[meta.vz[i_z],:]
        except:
            print 'Fatal error, possibly due to a too low turbulence intensity with respect to the demanded mean wind speed, try increase the input TI'
        DWM_TI_DATA = mfor.TI_DWM[meta.vz[i_z],:]
        ### Correct DWM_TI_DATA so that no point to have lower TI than "TIamb"
        DWM_TI_DATA[DWM_TI_DATA < np.nanmean(meta.mean_TI_DWM)] = np.nanmean(meta.mean_TI_DWM)
        for i_t in np.arange(0,len(meand.time),1):
            Ro_x                    = meand.meand_pos_x[i_z,i_t]
            Ro_y                    = meand.meand_pos_y[i_z,i_t]

            r_dist                  = np.sqrt((meta.x_mat - Ro_x)**2 + (meta.y_mat - Ro_y)**2 )
            tmp_index               = r_dist < mfor.WakeW[meta.vz[i_z]]*1.5
            tmp_field_WS            = np.ones((meta.nx,meta.ny))

            tmp_field_WS[tmp_index] = np.interp( r_dist[tmp_index],meta.vr_m, DWM_WS_DATA)
            ffor.WS_axial_ffor_tmp[:, :, i_z]  = ffor.WS_axial_ffor_tmp[:, :, i_z] + (tmp_field_WS)
            ffor.ffor_flow_field_ws_tmp2[:, :, i_z] = ffor.ffor_flow_field_ws_tmp2[:, :,i_z]+ (tmp_field_WS**2)

            tmp_field_TI            = meta.TI * np.ones((meta.nx,meta.ny))
            tmp_field_TI[tmp_index] = np.interp( r_dist[tmp_index],meta.vr_m,DWM_TI_DATA)

            ffor.ffor_flow_field_TI_tmp_tmp[:, :]      = tmp_field_TI
            ffor.TI_axial_ffor_tmp[:, :, i_z]     = ffor.TI_axial_ffor_tmp[:, :, i_z] + ffor.ffor_flow_field_TI_tmp_tmp**2

    # Stores the mean field
    for i_z in np.arange(0,meta.nz,1):
        ffor.TI_meand_axial_ffor[:, :, i_z]=np.sqrt(abs(ffor.ffor_flow_field_ws_tmp2[:, :, i_z] - ((ffor.WS_axial_ffor_tmp[:, :, i_z]**2)/len(meand.time)) )/ (len(meand.time)-1.0))
        ffor.WS_axial_ffor[:, :, i_z]      = (ffor.WS_axial_ffor_tmp[:, :, i_z]  / len(meand.time))
        ffor.TI_axial_ffor[:, :, i_z]      = (ffor.TI_axial_ffor_tmp[:, :, i_z]  / len(meand.time))**(1.0/2.0)

    # Store the ffor flow field
    ffor.x_vec                   = (meta.x_vec-meta.hub_x[0])/2.
    ffor.y_vec                   = (meta.y_vec-meta.hub_y)/2.
    ffor.z_vec                   = meta.z_vec+np.hstack((0., np.cumsum(meta.hub_z[0:])))[0]
    ffor.x_mat                   = meta.x_mat/2.
    ffor.y_mat                   = meta.y_mat/2.
    ffor.z_mat                   = meta.z_mat/meta.dz

    # print ffor.WS_axial_ffor.shape
    return ffor,meta

def DWM_make_grid(meta):
    """
    Function to adapt the grid (polar MFoR and cartesian FFoR) length
    in the streamwise direction based on the distance to the next turbine
    It also create extraction plane in the cartesian ffor domain

    Parameters
    ----------
    meta: class holding grid and ambient parameters

    Returns
    -------
    meta (updated)
    lz_mixl  : length of the mixing length domain
    vz_mixl  : vector of streamwise grid points of the polar MFoR domain in R
    vz       : vector of streamwise grid points indices of the cartesian domain
    nz       : vector length of streamwise grid points of the cartesian domain in D
    z_vec    : vector length of streamwise grid points of the cartesian domain in D
    """

    meta.vz=meta.hub_z[0:]*meta.dz
    # The mixL domain
    meta.lz_mixl=1.0+2.0*(max(meta.hub_z)) # in R: 1R longer than ffor flow field due to backward diff scheme
    meta.vz_mixl = np.linspace(0,meta.lz_mixl-meta.dz_mixl,meta.lz_mixl/meta.dz_mixl) # coordinate streamwise in R vector mixL

    meta.nz = len(meta.vz)           # nb points in z (streamwise) direction ffor flow field
    meta.z_vec= meta.vz/meta.dz

    return meta

def DWM_meta_meand(meand,meta):
    """
    Function to calculate the standard deviation of the wake center in-plane.

    Parameters
    ----------
    meand:   class holding the meandering parameters
    meta:    class holding grid and ambient parameters

    Returns
    -------
    meand: (updated)
    meand_pos_x [z,t]  : vector holding wake position in lateral
    meand_pos_y [z,t]  : vector holding wake position in longitunal
    std_meand_x [z,t]  : standard deviation of wake center position in lateral
    std_meand_y [z,t]  : standard deviation of wake center position in long
    time   [t]   : time vector
    """
    #meand.time             = np.arange(1,600.+1,1) # number of wake meandering samples 10 min at 1Hz
    meand.time             = np.arange(1,300+1,1) # number of wake meandering samples 10 min at 1Hz

    # Wake meandering based on meandering method with 0.8 as wake transport speed, see Keck et al. [4] & [5]
    meand.std_meand_x, meand.std_meand_y = meand_table_DWM_method(meta)

    # build meandering vectors at the planes specified by meta.vz
    meand.meand_pos_x=np.zeros((meta.nz,len(meand.time)))
    meand.meand_pos_y=np.zeros((meta.nz,len(meand.time)))

    # here the meandering patterns are assumed independent from each other between downstream planes
    # this is not physical as meandering paths are well correlated between consecutives downstream distances
    # with a scaling = variation of std_meand_x
    # time shift due to advection time

    seed_x=np.random.randn(len(meand.time),1)
    seed_y=np.random.randn(len(meand.time),1)


    for i_z in np.arange(0,meta.nz,1):
        meand.meand_pos_x[i_z,:]=(meta.hub_x[0] +  (meand.std_meand_x[i_z] *seed_x)).ravel()
        meand.meand_pos_y[i_z,:]=(meta.hub_y + (meand.std_meand_y[i_z] * seed_y)).ravel()
    return meand

def meand_table_DWM_method(meta):
    """
    Function to determine the meandering magnitude as function of ambient conditions at given downstream position i_z
    The parameters required for the determination of the meandering magnitude are: downstream distance, height,
    turbulence intensity and atmospheric stability
    
    Parameters
    ----------
    meta:    class holding grid and ambient parameters

    Returns
    -------
    std_meand_x[i_z,1]:   standard deviation of lateral meandering magnitude
    std_meand_y[i_z,1]:   standard deviation of longitudinal meandering magnitude

    """
    tmp_TI=np.zeros((1,2)).ravel();tmp_iTI=np.zeros((1,2)).ravel()
    location=meta.z_vec
    lloc=len(location)
    index_orig=np.argsort(location)
    location=np.sort(location) # sorting for interpolation on continuous function
    Meand = io.loadmat('../data/meand_data.mat')
    TI_vector     = np.array([0.0, 60.0, 100.0, 140.0, 200.0])/1000.0
    #creates tmp_TI for interpolation
    ind=TI_vector<=meta.TI
    r = np.array(range(len(ind)))
    tmp_iTI[0] = max(r[ind])
    tmp_TI[0]  = max(TI_vector[r[ind]])
    try:
        ind=TI_vector>=meta.TI
        r = np.array(range(len(ind)))
        tmp_iTI[1]= min(r[ind])
        tmp_TI[1]  = min(TI_vector[r[ind]])
    except:
        tmp_TI[1]  = tmp_TI[0]
    # Test if TI over TI table limit
    if (tmp_TI[0] == tmp_TI[1])==1:
        tmp_TI_length = 1
    else:
        tmp_TI_length = 2

    std_hor_tmp=np.zeros((tmp_TI_length,lloc));std_vert_tmp=np.zeros((tmp_TI_length,lloc))
    std_meand_x=np.zeros((lloc));std_meand_y=np.zeros((lloc))
    #std_hor_tmp=[];std_vert_tmp=[]
    for i_iTI in np.arange(0,tmp_TI_length,1):
        # Atmospheric stability effects based on Penas Mann stability & MARDM integration
        #According to 2013-01-11 mail
        tmp_data_hor= Meand[meta.atmo_stab][0][0][0][0][0][tmp_iTI[i_iTI]]
        tmp_data_vert= Meand[meta.atmo_stab][0][0][1][0][0][tmp_iTI[i_iTI]]

        #  Adds value at D=0 and 108D (=18D*6)
        tmp_data_hor  = np.vstack(([0.0, 0.0, 0.0],tmp_data_hor,tmp_data_hor[-1,:]*6.))
        tmp_data_vert = np.vstack(([0.0, 0.0, 0.0],tmp_data_vert,tmp_data_vert[-1,:]*6.))

        #specific the format of the meandering data
        dist_vector   = np.array([0.0, 1.0, 2.0, 3.0, 4.5, 6.0, 7.5, 9.0, 12.0, 18.0, 108.0])
        height_vector = np.array([40.0, 100.0, 160.0])

        #finds the wake meandering
        sp = interpolate.RectBivariateSpline(dist_vector, height_vector, tmp_data_hor, kx=1, ky=1, s=0)
        std_hor_tmp[i_iTI,:]=sp(location,meta.WTG_spec.H).reshape(1,lloc)

        sp = interpolate.RectBivariateSpline(dist_vector, height_vector, tmp_data_vert, kx=1, ky=1, s=0)
        std_vert_tmp[i_iTI,:]=sp(location,meta.WTG_spec.H).reshape(1,lloc)

    if (tmp_TI_length == 1)==1 and (tmp_TI[1] == meta.TI)==1: # use tmp_TI(1) value
        std_meand_x = std_hor_tmp[0,:]
        std_meand_y = std_vert_tmp[0,:]
    elif (tmp_TI_length == 1)==1: # Scale tmp_TI(1) value to TI_AMB value
        std_meand_x  = (meta.TI/tmp_TI[0])* std_hor_tmp[0,:]
        std_meand_y = (meta.TI/tmp_TI[0])* std_vert_tmp[0,:]
    else: # interpolate between tmp_TI(1) and tmp_TI(2).
        for i_z in np.arange(0,lloc,1):
            std_meand_x[i_z]  = np.interp(meta.TI,[tmp_TI[0], tmp_TI[1]],std_hor_tmp[:,i_z])
            std_meand_y[i_z] = np.interp(meta.TI,[tmp_TI[0], tmp_TI[1]],std_vert_tmp[:,i_z])

    # print std_meand_x
    std_meand_x=std_meand_x[index_orig] # reorder to original sorting
    std_meand_y=std_meand_y[index_orig]
    # print std_meand_x
    return std_meand_x, std_meand_y

def DWM_outputs(DWM,ffor,mfor,meta, aero, BEM):
    """
    Function that store chosen flow field and turbine results

    Parameters
    ----------
    DWM: list[nWT]
        list container for all results
    meta: class instance holding grid and ambient parameters
    ffor: class instance holding the fixed frame of reference flow field data
    mfor: class instance holding the meandering frame of reference flow field data

    Returns
    -------
    DWM: (updated)
       For the description of the returned variables, refer to the class definition.

    """
    dwm = Outputs()
    ###############  Store Flow field #################################################################################

    # Global flow field
    dwm.WS_axial_ffor          = ffor.WS_axial_ffor
    dwm.TI_axial_ffor          = ffor.TI_axial_ffor
    dwm.TI_meand_axial_ffor    = ffor.TI_meand_axial_ffor
    dwm.TI_tot_axial_ffor      = np.sqrt(ffor.TI_axial_ffor**2 + ffor.TI_meand_axial_ffor**2)
    # Store Mean turbine data
    dwm.WS_DWM                  =  meta.mean_WS_DWM
    dwm.TI_DWM                  =  meta.mean_TI_DWM
    dwm.C2C                     =  meta.C2C
    dwm.x_vec                   =  ffor.x_vec
    dwm.y_vec                   =  ffor.y_vec
    dwm.z_vec                   =  ffor.z_vec
    dwm.x_mat                   =  ffor.x_mat
    dwm.y_mat                   =  ffor.y_mat
    dwm.z_mat                   =  ffor.z_mat
    dwm.vr_mixl = meta.vr_mixl
    dwm.vz_mixl = meta.vz_mixl

    # Aero - BEM
    dwm.mean_a=aero.mean_a
    dwm.f_w=aero.f_w
    dwm.r_w=aero.r_w
    # Boundary conditions for U_w
    dwm.U_w=aero.U_w

    # Meandering frame of reference
    dwm.U_init= mfor.U_init
    dwm.Shear_add_du_dr=mfor.Shear_add_du_dr
    dwm.Shear_add_du_dz=mfor.Shear_add_du_dz
    dwm.TI_DWM=mfor.TI_DWM
    dwm.U=mfor.U
    dwm.U_init=mfor.U_init
    dwm.V=mfor.V
    dwm.WakeW=mfor.WakeW
    dwm.du_dr_DWM=mfor.du_dr_DWM
    dwm.du_dr_tot=mfor.du_dr_tot
    dwm.visc=mfor.visc

    ################### Save to DWM list ###############################################################
    DWM[str(meta.wtg_ind[0])]=dwm

    return DWM
