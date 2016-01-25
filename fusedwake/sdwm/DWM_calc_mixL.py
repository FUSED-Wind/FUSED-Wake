# -*- coding: utf-8 -*-
""" Python interpretation of Rolf-Erik Keck's Ainslie model
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
from math import pi
from scipy import linalg

def DWM_init_calc_mixl(meta,mfor):
    """Function that initializes the Ainslie model including improvements from
    [Keck_2011]_, [Keck_2012]_, [Keck_2013a]_, [Keck_2013b]_, [Keck_2014]_,
    [Keck_2015]_

    Parameters
    ----------
    meta : (instance of class)
        Instance of class Meta holding DWM core variables
    mfor : (instance of class)
        Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model

    Returns
    -------
    mfor : (instance of class)
        Update instance of class Mfor with inialized flow variables based on wake case
    F1_vector : (np.array float)
        DWM filter functions F 1 governing the development of turbulent stresses [Keck_2011]_, [Keck_2012]_, eq 3 of [Keck_2015]_
    F2_vector : (np.array float)
        DWM filter functions F 2 governing the development of turbulent stresses [Keck_2011]_, [Keck_2012]_, eq 3 of [Keck_2015]_
    visc_wake1 : (np.array float)
        initialized contribution due to ambient turbulence for the eddy viscosity [Keck_2015]_, Eq 3 first term of right and side
    visc_wake2 : (np.array float)
        initialized contribution due to shear layer of wake deficit for the eddy viscosity [Keck_2015]_, Eq 3 last term of right and side
    visc_wake : (np.array float)
        eddy viscosity of DWM, combination of visc_wake1 and visc_wake2
    u_star_DEF : (np.array float)
        non dimensionalivelocity scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
    l_star_DEF : (np.array float)
        non dimensional integral length scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
    One_div_du_dr_DWW : (np.array float)
        denominator of Eq (5) in [Keck_2015]_
    width : (np.array int)
        wake width vector index
    """
    # Generate the DWM filter functions for the eddy viscosity formulation
    F1_vector  = np.hstack((np.linspace(meta.f1[0],1,meta.f1[1]*meta.dz/2),np.linspace(1,1,meta.lz_mixl*meta.dz/2)))
    F2_z_vec   = np.arange(2+1./meta.dz,(len(F1_vector)+1)*(1./meta.dz),1./meta.dz)
    F2_vector  = np.hstack((np.linspace(meta.f2[0],meta.f2[0],2*meta.dz), 1.-(1.-meta.f2[0])*np.exp(-meta.f2[1]*(F2_z_vec-2.))))

    # initiate the U, V, visc etc... matrices
    mfor.V                 = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float)
    mfor.U                 = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.visc              = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.du_dr_DWM         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.du_dr_tot         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.Turb_Stress_DWM   = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.TI_DWM            = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake1        = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake2        = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    One_div_du_dr_DWM = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )

    # DWM boundary conditions
    # Rotor plane
    mfor.U[0,:] = mfor.U_init
    #Centerline
    mfor.V[0,:]     = 0.0
    # Atmospheric stability effects based on Keck et al. [Keck_2013a]_
    if (meta.atmo_stab== 'VU')==1:
          L_ABL_vector         = [42.9157,      68.5912,      88.0709]
          UW_UU_vector         = [-0.27991,    -0.26012,    -0.23296]
          L_DEF_vector         = [14.8641,      19.5041,      22.0656]
          UU_DEF_UU_ABL_vector = [0.51678,     0.47861,     0.42581]
          UW_DEF_UU_DEF_vector = [-0.23688,    -0.13927,   -0.097113]
    elif (meta.atmo_stab== 'U')==1:
          L_ABL_vector         = [39.0648,      61.8843,      80.8478]
          UW_UU_vector         = [-0.27441,    -0.28051,    -0.26325]
          L_DEF_vector         = [14.0137,      18.1673,      20.8881]
          UU_DEF_UU_ABL_vector = [0.51067,      0.4433,     0.41874]
          UW_DEF_UU_DEF_vector = [-0.2544,    -0.18689,    -0.12968]
    elif (meta.atmo_stab== 'NU')==1:
          L_ABL_vector         = [33.2038,       48.321,      66.7102]
          UW_UU_vector         = [-0.27136,    -0.28125,    -0.28078]
          L_DEF_vector         = [12.7207,      15.9199,      18.7976]
          UU_DEF_UU_ABL_vector = [0.54998,     0.49571,      0.4139]
          UW_DEF_UU_DEF_vector = [-0.26641,    -0.22075,    -0.18199]
    elif (meta.atmo_stab== 'N')==1:
          L_ABL_vector         = [26.5352,       34.026,      40.7458]
          UW_UU_vector         = [-0.27359,    -0.27887,    -0.27935]
          L_DEF_vector         = [11.065,      12.9746,      14.4395]
          UU_DEF_UU_ABL_vector = [0.63044,     0.57982,      0.5287]
          UW_DEF_UU_DEF_vector = [-0.27341,    -0.25684,    -0.24217]
    elif (meta.atmo_stab== 'NS')==1:
          L_ABL_vector         = [21.2064,      27.1416,      34.2689]
          UW_UU_vector         = [-0.27331,    -0.27636,    -0.28028]
          L_DEF_vector         = [9.50836,      11.2453,      13.0561]
          UU_DEF_UU_ABL_vector = [0.69202,     0.63823,     0.59067]
          UW_DEF_UU_DEF_vector = [-0.27954,    -0.26985,    -0.25258]
    elif (meta.atmo_stab== 'S')==1:
          L_ABL_vector         = [12.4648,      14.6943,      21.5762]
          UW_UU_vector         = [-0.27085,    -0.27464,     -0.2781]
          L_DEF_vector         = [6.3775,      7.2553,      9.6425]
          UU_DEF_UU_ABL_vector = [0.80217,     0.78088,     0.71086]
          UW_DEF_UU_DEF_vector = [-0.28288,     -0.2835,     -0.2758]
    elif (meta.atmo_stab== 'VS')==1:
          L_ABL_vector         = [7.4849,      10.4295,      23.5966]
          UW_UU_vector         = [-0.26527,     -0.2729,    -0.26608]
          L_DEF_vector         = [4.21111,      5.53777,      10.2117]
          UU_DEF_UU_ABL_vector = [0.87381,     0.83775,     0.63848]
          UW_DEF_UU_DEF_vector = [-0.27585,    -0.28304,    -0.27871]
    L_ABL          = np.interp(meta.WTG_spec.H,[40.,  100.,  160.], L_ABL_vector  )    # int_Lww(k3) i.e (Integral length scale (in vertical directions), from ww(k3))
    UW_UU          = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UW_UU_vector )    # ratio of UW and UU stresses for whole spectra
    L_DEF          = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], L_DEF_vector )   #    Integral length scale (in vertical directions), Meandering length scale subtracted, from ww(k3)
    UU_DEF_UU_ABL  = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UU_DEF_UU_ABL_vector)    # Part of normal stress in the deficit module
    UW_DEF_UU_DEF  = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UW_DEF_UU_DEF_vector)    # ratio of UW and UU stres

    Rotor_R        = 40.   # DO NOT CHANGE!!  #ATMOSTAB ANALYSIS IS CARRIED OUT OVER R = 40m, which should be used to normalize the length scales
    l_star_ABL     = L_ABL / Rotor_R#
    l_star_DEF     = L_DEF / Rotor_R#

    # Normalize UU_160m to neutral condition (to use calibration from Keck et al. [3])
    UU_DEF_UU_ABL_fac  = np.interp(meta.WTG_spec.H,[40., 100., 160.], [0.63044,     0.57982,      0.5287])
    UU_DEF_UU_ABL      = UU_DEF_UU_ABL / UU_DEF_UU_ABL_fac

    #CALCULATE u* according to:
    #1. u* ~= (mean(u'w')^2 )^0.25
    #2. {mean(u'w') = mean(u'u')*Cuw_uu}
    #3. {u' ~= TI (in normalized form)}
    # => u* ~= ((TI^2 * Cuw_uu )^2)^0.25
    u_star_ABL     = ( (  (meta.TI)**2            * np.abs(UW_UU)         )**2 )**0.25 # replaced by inflow TI
    u_star_DEF     = ( (  (meta.mean_TI_DWM)**2 * UU_DEF_UU_ABL * np.abs(UW_DEF_UU_DEF) )**2 )**0.25
    mfor.Shear_add_du_dz = u_star_ABL / l_star_ABL
    width=np.zeros((len(np.arange(1,len(meta.vz_mixl),1)), 1),dtype=int)
    return mfor,F1_vector,F2_vector,visc_wake1,visc_wake2,visc_wake,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,width


def DWM_calc_wake_width(mfor,width,meta,j):
    """Function that estimates the wake width.

        Parameters
        ----------
        meta : (instance of class)
            Instance of class Meta holding DWM core variables
        mfor : (instance of class)
            Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        width (np.array int) : wake width vector index
        j (int) : index in the main Ainslie forward scheme loop

        Returns
        -------
        width : (np.array int)
            updated wake width vector index
    """

    # Calculating wake width
    dr_DWM        = 1.0/meta.dr_mixl
    r_vec_DWM     = np.arange(dr_DWM/2.0,meta.lr_mixl-dr_DWM/2.0+dr_DWM,dr_DWM)
    r_vec_DWM = r_vec_DWM[r_vec_DWM<meta.lr_mixl] # bug fix
    dA_DWM        = np.array([pi*r_vec_DWM[1:]**2-pi*r_vec_DWM[0:-1]**2])
    b_loop=meta.dr_mixl

    np.linspace(0,meta.lr_mixl-meta.dr,(meta.lr_mixl)/meta.R_WTG*meta.dr_mixl).shape
    Def_DWM       = np.sum(( 1.0-mfor.U[j-1,1:])*dA_DWM  )

    while(b_loop):
        b_loop = b_loop +1
        Def_DWM_mixL = np.sum((1.0 - mfor.U[j-1,1:(b_loop)]) * dA_DWM[0,0:(b_loop-1)])
        if (Def_DWM_mixL > Def_DWM * 0.95)==1:
            break
        elif (b_loop == (meta.dr_mixl * meta.lr_mixl) -1)==1:
            break

    width[j-1]    =  b_loop
    return width

def DWM_eddy_viscosity(mfor,meta,width,visc_wake,visc_wake1,visc_wake2,F1_vector,F2_vector,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,j):
    """Function that calculates the eddy viscosity

    Parameters
    ----------
    meta : (instance of class)
        Instance of class Meta holding DWM core variables
    mfor : (instance of class)
        Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
    F1_vector : (np.array float)
        DWM filter functions F 1 governing the development of turbulent stresses [Keck_2011]_, [Keck_2012]_, eq 3 of [Keck_2015]_
    F2_vector : (np.array float)
        DWM filter functions F 2 governing the development of turbulent stresses [Keck_2011]_, [Keck_2012]_, eq 3 of [Keck_2015]_
    visc_wake1 : (np.array float)
        initialized contribution due to ambient turbulence for the eddy viscosity [Keck_2015]_, Eq 3 first term of right and side
    visc_wake2 : (np.array float)
        initialized contribution due to shear layer of wake deficit for the eddy viscosity [Keck_2015]_, Eq 3 last term of right and side
    visc_wake : (np.array float)
        eddy viscosity of DWM, combination of visc_wake1 and visc_wake2
    u_star_DEF : (np.array float)
        non dimensional velocity scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
    l_star_DEF : (np.array float)
        non dimensional integral length scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
    One_div_du_dr_DWW : (np.array float)
        denominator of Eq (5) in [Keck_2015]_
    width : (np.array int)
        wake width vector index
    j (int) : index in the main Ainslie forward scheme loop

    Returns
    -------
    mfor : (instance of class)
        updated instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
    """

## Calculate eddy viscosity
    ## Include blend between original Prandtl model and Ainslie to avoid issues when wake turbulence goes to 0.
    ## The largest eddy viscosity at each point is applied.
    # Calculate mean flow gradient - du/dr is created with CDS (apart from 1st and last point)
    mfor.du_dr_DWM[j-1,0]                  = (mfor.U[j-1,1] - mfor.U[j-1,0])/meta.dr
    mfor.du_dr_DWM[j-1,1:meta.lr_mixl*meta.dr_mixl-2]  = (mfor.U[j-1,2:(meta.lr_mixl*meta.dr_mixl-2)+1] \
    - mfor.U[j-1,0:(meta.lr_mixl*meta.dr_mixl-2)-1])/(2*meta.dr)
    mfor.du_dr_DWM[j-1,meta.lr_mixl*meta.dr_mixl-1]      = (mfor.U[j-1,meta.lr_mixl*meta.dr_mixl-1] \
    - mfor.U[j-1,meta.lr_mixl*meta.dr_mixl-2])/meta.dr

    # Blend of mixL and Ainslie eddy visc
    visc_wake1[j-1,:]     = F2_vector[j-1]* meta.k2 *( meta.vr_mixl[width[j-1]-1]/meta.R_WTG )**2 * np.abs(mfor.du_dr_DWM[j-1,:])
    visc_wake2[j-1,:]     = F2_vector[j-1]* meta.k2 *( meta.vr_mixl[width[j-1]-1]/meta.R_WTG )   * (1.0 - np.min(mfor.U[j-1,:]) )
    visc_wake[j-1,:]      = np.maximum(visc_wake1[j-1,:],visc_wake2[j-1,:])
    #max operator is included in the eddy viscosity formulation to avoid underestimating the
    #turbulent stresses at locations where the velocity gradient of the deficit du_dr approaches zero
    # Atmospheric eddy visc as u*l*, yields total eddy viscosity
    visc_norm_factor      = 6.3918 # Applied to use Keck et al. [3] calibration
    mfor.visc[j-1,:]           = F1_vector[j-1]*meta.k1*visc_norm_factor*u_star_DEF*l_star_DEF + visc_wake[j-1,:]

    ## Include contribution from atmospheric boundary layer on DWM
    ##  turbulent stresses. This effect is taken into account by:
    # 1. Calculate the azimuthally averaged local gradient (du/dr tot) acting of the eddy viscosity as a combination of du/dr in the DWM model and du/dz from ABL
    # 2. The du/dr contribution is constant in azimuthal direction. The du/dz part is assumed linear, which gives a sinus curve in a du/dr system
    # 3. Instead of manipulating the velocity field, the eddy viscosity is increased by a "du/dr_total / du/dr_DWM"
    #=> Visc*        =  Visc * du/dr_total / du/dr_DWM
    #   => Turb_Stress  =  Visc* * du/dr_DWM = Visc * (du/dr_total / du/dr_DWM) * du/dr_DWM =  Visc * du/dr_total
    # 4. "Wiener filter" is used to avoid problems when du/dr = 0, idea:  1/f(x) ~= f(x) / (f(x)^2 + k)

    # Calculate total mean flow gradient - adds shear contribution via
    # sinus function. This gets the stresses right, but sign is wrong in
    #regions where du/dr_DWM - sign of du/dz_ABL is negative
    # notations as per Keck et al. [3].
    du_dr_DWM_du_dz=np.array(np.absolute(mfor.du_dr_DWM[j-1,:]) / mfor.Shear_add_du_dz ,dtype=complex)
    alfa_1      = np.arcsin(du_dr_DWM_du_dz)
    alfa_2      = pi - alfa_1
    alfa_1=np.asarray([abs(x) for x in alfa_1])
    alfa_2=np.asarray([abs(x) for x in alfa_2])
    mfor.du_dr_tot[j-1,0:meta.lr_mixl*meta.dr_mixl] = ( np.absolute(mfor.du_dr_DWM[j-1,:]) *2.0*pi +\
    ((np.absolute(mfor.du_dr_DWM[j-1,:]) < mfor.Shear_add_du_dz) * 2.0 * \
    (mfor.Shear_add_du_dz*2.0*np.cos(alfa_1) - np.absolute(mfor.du_dr_DWM[j-1,:])*(alfa_2 - alfa_1) ) ) ) / (2.0*pi)
    # du/dr_DWM block of area
    # condition for added shear gradient (if du/dr_DWM >= du/dz_ABL there are no contribution)
    # Area A1 + A2 in figure XXX
    # Scaling from area to gradient
    k_wiener                 = 2.0*mfor.Shear_add_du_dz * meta.dr**2
    One_div_du_dr_DWM[j-1,:] = mfor.du_dr_DWM[j-1,:] / (mfor.du_dr_DWM[j-1,:]**2 + k_wiener)
    visc_fac                 = np.maximum(1.0, (mfor.du_dr_tot[j-1,:] * np.fabs(One_div_du_dr_DWM[j-1,:])))
    mfor.visc[j-1,:]              = mfor.visc[j-1,:] * visc_fac

    return mfor


def DWM_velocity_solver(mfor,meta,j):
    """Function that defines tridiagonal matrix to solve the NS equations Eq 1 and in [Keck_2015]_
       The momentum equation is discretized using a second order central difference scheme in radial direction and a first order upwind scheme in flow direction

    (1) Solve the momentum equation for the streamwise velocity component at all radial positions explicitly, by using the value of the radial velocity component and the eddy viscosity from the previous location upstream. This yields a tri- diagonal equation system where all the coefficients are known, which can easily be solved by any tridiagonal ma- trix algorithm.
    (2) Once the streamwise velocity is known, the radial velocity for all radial positions can be updated using the continuity equation
    (3) The eddy viscosity for all radial positions is updated using Eq. (3) in [Keck_2011]_
    (4) March to the next downstream location and repeat steps 1-3.

    Parameters
    ----------
    meta : (instance of class)
        Instance of class Meta holding DWM core variables
    mfor : (instance of class)
        Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
    j (int) : index in the main Ainslie forward scheme loop

    Returns
    -------
    HL (np.array float)
    mat (np.array float)

    """

    # init radial vector
    ind_R=range(1,int(meta.lr_mixl*meta.dr_mixl-1),1)
    ind_R_p=[z+1 for z in ind_R]
    ind_R_m=[z-1 for z in ind_R]

    HL = np.zeros((meta.lr_mixl*meta.dr_mixl))
    mat = np.zeros(((meta.lr_mixl*meta.dr_mixl),(meta.lr_mixl*meta.dr_mixl)))

    # Input BC for wake center
    HL[0]     = (mfor.U[j-1,0]**2     / meta.dz_mixl)
    HL[ind_R] = (mfor.U[j-1,ind_R]**2/meta.dz_mixl)
    HL[meta.lr_mixl*meta.dr_mixl-1]       = (mfor.U[j-1,meta.lr_mixl*meta.dr_mixl-1]/ meta.dz_mixl)

    mat[0,0]  =  mfor.U[j-1,0]/meta.dz_mixl     + (2.0*mfor.visc[j-1,0]/(meta.dr**2))
    mat[1,0]  = -(2.0*mfor.visc[j-1,0] /(meta.dr**2))

    VL11=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL21=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL31=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL41=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL12=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL13=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL22=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL23=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL33=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL43=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL1=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL2=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)
    VL3=np.zeros((meta.lr_mixl*meta.dr_mixl)-1)

    # Calculation of U for the wake body
    # mfor.visc[j-1,range(0,int(meta.lr_mixl*meta.dr_mixl),1)+1]
    VL11[ind_R]             = -mfor.V[j-1,ind_R]      / (2.0*meta.dr)
    VL21[ind_R]             = mfor.visc[j-1,ind_R]    / (2.0*meta.vr_mixl[ind_R]*meta.dr)
    VL31[ind_R]             = -mfor.visc[j-1,ind_R]   / (meta.dr**2)
    VL41[ind_R]             = (mfor.visc[j-1,ind_R_p] - mfor.visc[j-1,ind_R_m])  / (2*meta.dr)**2 # new term due to d(nu_t)/dr dependence
    VL12[ind_R]             = mfor.U[j-1,ind_R]       / (meta.dz_mixl)
    VL22[ind_R]             = +2.0*mfor.visc[j-1,ind_R] / (meta.dr**2)
    VL13[ind_R]             = mfor.V[j-1,ind_R]       / (2.0*meta.dr)
    VL23[ind_R]             = -mfor.visc[j-1,ind_R]   / (2.0*meta.vr_mixl[ind_R]*meta.dr)
    VL33[ind_R]             = -mfor.visc[j-1,ind_R]   / (meta.dr**2)
    VL43[ind_R]             = -(mfor.visc[j-1,ind_R_p] - mfor.visc[j-1,ind_R_m])  / (2.0*meta.dr)**2 # new term due to d(nu_t)/dr dependence
    VL1[ind_R]              = VL11[ind_R] + VL21[ind_R] + VL31[ind_R] + VL41[ind_R]
    VL2[ind_R]              = VL12[ind_R] + VL22[ind_R]
    VL3[ind_R]              = VL13[ind_R] + VL23[ind_R] + VL33[ind_R] + VL43[ind_R]

    # build the matrix for X =A/B
    mat[ind_R_m,ind_R] = VL1[ind_R]
    mat[ind_R ,ind_R] = VL2[ind_R]
    mat[ind_R_p,ind_R] = VL3[ind_R]

    # Input BC for wake edge
    VL1                   = 0.0
    VL2                   = 1.0/meta.dz_mixl
    mat[meta.lr_mixl*meta.dr_mixl-2,meta.lr_mixl*meta.dr_mixl-1]  = VL1
    mat[meta.lr_mixl*meta.dr_mixl-1,meta.lr_mixl*meta.dr_mixl-1]  = VL2

    mat=mat.T
    HL=HL.T

    return HL, mat

def DWM_calc_mixL(meta, aero, mfor):
    """ Main Ainslie - mixing length [Keck_2011]_ function that computes the wake deficit as function of downstream distance in the meandering frame of reference

    (1) Solve the momentum equation for the streamwise velocity component at all radial positions explicitly, by using the value of the radial velocity component and the eddy viscosity from the previous location upstream. This yields a tri- diagonal equation system where all the coefficients are known, which can easily be solved by any tridiagonal matrix algorithm.
    (2) Once the streamwise velocity is known, the radial velocity for all radial positions can be updated using the continuity equation
    (3) The eddy viscosity for all radial positions is updated using Eq. (3) in [Keck_2011]_
    (4) March to the next downstream location and repeat steps 1-3

    Parameters
    ----------
        meta : (instance of class)
            Instance of class Meta holding DWM core variables
        aero : (instance of class)
            Instance of class Aero holding BEM-aero core variables
        mfor : (instance of class)
            Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model

    Returns
    -------
        mfor : (instance of class)
            Updated instance of class Mfor holding the velocity deficit

    """

    # Init arrays and constants for mfor_mixL module
    mfor,F1_vector,F2_vector,visc_wake1,visc_wake2,visc_wake,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,width=DWM_init_calc_mixl(meta,mfor)

    #  Start x-stepping & solving of BLE equations
    for j in np.arange(1,len(meta.vz_mixl),1):
        ## Calculating wake width
        width=DWM_calc_wake_width(mfor,width,meta,j)

        ## Calculate eddy viscosity
        mfor=DWM_eddy_viscosity(mfor,meta,width,visc_wake,visc_wake1,visc_wake2,F1_vector,F2_vector,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,j)

        ## Calculation matrix
        HL, mat=DWM_velocity_solver(mfor,meta,j)

        ## Solve for U
        mfor.U[j,:] = linalg.solve(mat, HL)

        for i in np.arange(0,meta.lr_mixl*meta.dr_mixl-1,1):
            mfor.V[j,i+1] = (meta.vr_mixl[i] / meta.vr_mixl[i+1]) * mfor.V[j,i] - (meta.dr/(2*meta.dz_mixl))*( (mfor.U[j,i+1] - mfor.U[j-1,i+1]) + \
            (meta.vr_mixl[i] / meta.vr_mixl[i+1])*((mfor.U[j,i] - mfor.U[j-1,i])) )

        # POST PROCESSING SIGNAL: Turbulent stress
        mfor.Turb_Stress_DWM[j-1,:]           = mfor.visc[j-1,:]  * mfor.du_dr_DWM[j-1,:]

        # POST PROCESSING SIGNAL: TI_DWM mforulated based on the relation derived between u'v' and u'u' (see Keck et al. [3])
        x_uw_wake      = 1.0
        C_uw_wake      = (0.7550 - meta.mean_TI_DWM*1.75) / 2 # Article states: "C_uw_wake = 0.3", but this seems to give better results (=> 0.3 for TI=8.6%)
        mfor.TI_DWM[j-1,:]  = np.sqrt( np.abs( (1.0 / (x_uw_wake * C_uw_wake)) * mfor.Turb_Stress_DWM[j-1,:]) )
        mfor.WakeW = meta.vr_mixl[width-1]

    return mfor
