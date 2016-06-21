"""GC Larsen wake model applied to offshore wind farms (WindFarm object)

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

References:
[1] Larsen GC. "A simple stationary semi-analytical wake model", 2009

"""
import numpy as np
import fusedwake.WindTurbine as wt
import fusedwake.WindFarm as wf
from copy import copy

def get_r96(D, CT, TI, pars=[0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0]):
    """Computes the wake radius at 9.6D downstream location of a turbine

    .. math::
        R_{9.6D} = a_1 \\exp (a_2 C_T^2 + a_3 C_T + a_4)  (b_1  TI + b_2)  D

    Inputs
    ----------
    D: float
        Wind turbine diameter
    CT: float
        Outputs WindTurbine object's thrust coefficient
    TI: float
        Ambient turbulence intensity
    pars: list
        GCL Model parameters [a1, a2, a3, a4, b1, b2]

    Returns
    -------
    R96: float
        Wake radius at 9.6D downstream location
    """
    a1, a2, a3, a4, b1, b2 = pars
    R96 = a1 * (np.exp(a2 * CT * CT + a3 * CT + a4)) * (b1 * TI + b2) * D

    return R96

def get_Rw(x, R, TI, CT, pars=[0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0]):
    """Computes the wake radius at a location.
    [1]-eq.3

    .. math::
        R_w = \\left(\\frac{105  c_1^2 }{2 \\pi}\\right)^{0.2} (C_T A (x + x_0))^{1/3}

    with A, the area, and x_0 and c_1 defined as

    .. math::
        x_0 = \\frac{9.6 D}{\\left(\\frac{2 R_96}{k D} \\right)^3 - 1}

        c_1 = \\left(\\frac{k D}{2}\\right)^{5/2}
              \\left(\\frac{105}{2 \\pi} \\right)^{-1/2}
              (C_T A x_0)^{-5/6}

    with k and m defined as

    .. math::
        k = \\sqrt{\\frac{m + 1}{2}}

        m = \\frac{1}{\\sqrt{1 - C_T}}

    Inputs
    ----------
    x: float or ndarray
        Distance between turbines and wake location in the wind direction
    R: float
        Wind turbine radius
    TI: float
        Ambient turbulence intensity
    CT: float
        Outputs WindTurbine object's thrust coefficient

    Returns
    -------
    Rw: float or ndarray
        Wake radius at a location
    """
    _ones = np.ones(np.shape(x))
    D = 2.0 * R
    Area = np.pi * D**2.0 / 4.0

    m = 1.0 / (np.sqrt(1.0 - CT))
    k = np.sqrt((m + 1.0) / 2.0)

    R96 = get_r96(D, CT, TI, pars)

    x0 = (9.6 * D) / ((2.0 * R96 / (k * D))**3.0 - 1.0)
    term1 = (k * D / 2.0)**2.5
    term2 = (105.0/(2.0*np.pi))**-0.5
    term3 = (CT * Area * x0)**(-5.0 / 6.0)
    c1 = term1 * term2 * term3

    Rw = ((105.0 * c1**2.0 / (2.0 * np.pi))**0.2) * (CT * Area * (x + x0 * _ones))**(1.0 / 3.0)

    if type(x) == float and x+x0 <= 0.: Rw = 0
    elif type(x) == np.ndarray: Rw[x + x0 * _ones <= 0.] = 0.
    return Rw, x0, c1

def get_dU(x,r,R,CT,TI,
    order=1,
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake velocity deficit at a location

    Inputs
    ----------
    x: float
        Distance between turbines and wake location in the wind direction
    r: float
        Radial distance between the turbine and the location
    R: float
        Wake producing turbine's radius [m]
    CT: float
        Outputs WindTurbine object's thrust coefficient
    TI: float
        Ambient turbulence intensity [-]
    order: int, optional

    Returns
    -------
    dU: float
        Wake velocity deficit at a location
    """
    _ones = np.ones(np.shape(x))

    D = 2.*R
    Area=np.pi*D*D/4.
    Rw, x0, c1 = get_Rw(x, R, TI, CT, pars)

    term10=(1./9.)*_ones
    term20=(CT*Area*(x+x0)**(-2.))**(1./3.)
    term310=(r**(1.5))
    term320=(3.*c1*c1*CT*Area*(x+x0))**(-0.5)
    term30=term310*term320
    term41=(35./(2.*np.pi))**(3./10.)
    term42=(3.*c1*c1)**(-0.2)
    term40=term41*term42
    dU1=-term10*term20*(term30-term40)**2.

    dU = dU1

    if order == 2:

        z_term1 = r**1.5
        z_term2 = (CT*Area*(x+x0))**(-0.5)
        z_term3 = ((35./(2.*np.pi))**(-3./10.))*((3.*c1*c1)**(-3./10.))
        z = z_term1*z_term2*z_term3

        d_term = (4./81.)*(((35./(2.*np.pi))**(6./5.))*((3.*c1*c1)**(-12./15.)))

        d_4_const = (1./40.)
        d_3_const = (-4.+48./40.)*1./19.
        d_2_const = (6.+27.*d_3_const)*1./4.
        d_1_const = (4.-12.*d_2_const)*1./5.
        d_0_const = (-1.-3.*d_2_const)*1./8.

        d_0 = d_term*d_0_const
        d_1 = d_term*d_1_const
        d_2 = d_term*d_2_const
        d_3 = d_term*d_3_const
        d_4 = d_term*d_4_const

        dU2_const = ((CT*Area*((x+x0)**(-2.)))**(2./3.))
        dU2_term0 = d_0*(z**0.)
        dU2_term1 = d_1*(z**1.)
        dU2_term2 = d_2*(z**2.)
        dU2_term3 = d_3*(z**3.)
        dU2_term4 = d_4*(z**4.)

        dU2 = dU2_const*(dU2_term0+dU2_term1+dU2_term2+dU2_term3+dU2_term4)

        dU=dU1 + dU2

    if type(r)==np.ndarray: dU[Rw<r]=0. # Outside the wake
    elif type(r)==float and Rw<r: dU = 0.

    if type(x)==np.ndarray: dU[x<=0.]=0. # upstream the wake gen. WT
    elif type(x)==float and x<=0.: dU = 0.

    if CT==0: dU = 0.0*dU

    return dU

def get_dUeq(x,y,z,RT,R,CT,TI,
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake velocity deficit at a location

    Inputs
    ----------
    x: float or array
        Distance between wake generating turbines and wake operating
        turbines in the streamwise direction
    y: float or array
        Distance between wake generating turbines and wake operating
        turbine in the crossflow horizontal direction
    z: float or array
        Distance between wake generating turbines and wake operating
        turbine in the crossflow vertical direction
    RT: float
        Wake operating turbine's radius [m]
    R: float or array
        Wake generating turbine's radius [m]
    TI: float
        Ambient turbulence intensity for the wake generating turbine [-]
    CT: float
        Thrust coefficient for the wake generating turbine [-]
    order: int, optional

    Returns
    -------
    dUeq: float
        Rotor averaged wake velocity deficit for each wake operating WT
    """

    # New improved quadrature rule for wake deficit rotor averaging
    node_R, node_th, weight = np.array([[ 0.26349922998554242692 ,  4.79436403870179805864 ,  0.00579798753740115753 ],
    [ 0.26349922998554242692 ,  5.13630491629471475079 ,  0.01299684397858970851 ],
    [ 0.26349922998554242692 ,  5.71955352542765460555 ,  0.01905256317618122044 ],
    [ 0.26349922998554242692 ,  0.20924454049880022999 ,  0.02341643323656225281 ],
    [ 0.26349922998554242692 ,  1.10309379714216659885 ,  0.02569988335562909190 ],
    [ 0.26349922998554242692 ,  2.03849885644762496284 ,  0.02569988335562912660 ],
    [ 0.26349922998554242692 ,  2.93234811309099407950 ,  0.02341643323656214179 ],
    [ 0.26349922998554242692 ,  3.70522443534172518653 ,  0.01905256317618119616 ],
    [ 0.26349922998554242692 ,  4.28847304447466459720 ,  0.01299684397858971198 ],
    [ 0.26349922998554242692 ,  4.63041392206758217753 ,  0.00579798753740114539 ],
    [ 0.57446451431535072718 ,  4.79436403870179805864 ,  0.01086984853977092380 ],
    [ 0.57446451431535072718 ,  5.13630491629471475079 ,  0.02436599330905551281 ],
    [ 0.57446451431535072718 ,  5.71955352542765460555 ,  0.03571902745281423097 ],
    [ 0.57446451431535072718 ,  0.20924454049880022999 ,  0.04390024659093685194 ],
    [ 0.57446451431535072718 ,  1.10309379714216659885 ,  0.04818117282305908744 ],
    [ 0.57446451431535072718 ,  2.03849885644762496284 ,  0.04818117282305915683 ],
    [ 0.57446451431535072718 ,  2.93234811309099407950 ,  0.04390024659093664378 ],
    [ 0.57446451431535072718 ,  3.70522443534172518653 ,  0.03571902745281418240 ],
    [ 0.57446451431535072718 ,  4.28847304447466459720 ,  0.02436599330905552321 ],
    [ 0.57446451431535072718 ,  4.63041392206758217753 ,  0.01086984853977089951 ],
    [ 0.81852948743000586429 ,  4.79436403870179805864 ,  0.01086984853977090992 ],
    [ 0.81852948743000586429 ,  5.13630491629471475079 ,  0.02436599330905548505 ],
    [ 0.81852948743000586429 ,  5.71955352542765460555 ,  0.03571902745281418934 ],
    [ 0.81852948743000586429 ,  0.20924454049880022999 ,  0.04390024659093679643 ],
    [ 0.81852948743000586429 ,  1.10309379714216659885 ,  0.04818117282305903193 ],
    [ 0.81852948743000586429 ,  2.03849885644762496284 ,  0.04818117282305909438 ],
    [ 0.81852948743000586429 ,  2.93234811309099407950 ,  0.04390024659093658826 ],
    [ 0.81852948743000586429 ,  3.70522443534172518653 ,  0.03571902745281413383 ],
    [ 0.81852948743000586429 ,  4.28847304447466459720 ,  0.02436599330905549199 ],
    [ 0.81852948743000586429 ,  4.63041392206758217753 ,  0.01086984853977088737 ],
    [ 0.96465960618086743494 ,  4.79436403870179805864 ,  0.00579798753740116100 ],
    [ 0.96465960618086743494 ,  5.13630491629471475079 ,  0.01299684397858971545 ],
    [ 0.96465960618086743494 ,  5.71955352542765460555 ,  0.01905256317618123432 ],
    [ 0.96465960618086743494 ,  0.20924454049880022999 ,  0.02341643323656226669 ],
    [ 0.96465960618086743494 ,  1.10309379714216659885 ,  0.02569988335562910925 ],
    [ 0.96465960618086743494 ,  2.03849885644762496284 ,  0.02569988335562914394 ],
    [ 0.96465960618086743494 ,  2.93234811309099407950 ,  0.02341643323656215567 ],
    [ 0.96465960618086743494 ,  3.70522443534172518653 ,  0.01905256317618120657 ],
    [ 0.96465960618086743494 ,  4.28847304447466459720 ,  0.01299684397858972065 ],
    [ 0.96465960618086743494 ,  4.63041392206758217753 ,  0.00579798753740114886 ]]).T

    x_msh, node_R_msh = np.meshgrid(x,node_R)
    y_msh, node_th_msh = np.meshgrid(y,node_th)
    z_msh, weight_msh = np.meshgrid(z,weight)

    xe = x_msh
    ye = y_msh + RT*node_R_msh*np.cos(node_th_msh)
    ze = z_msh + RT*node_R_msh*np.sin(node_th_msh)
    re = np.sqrt( ye**2. + ze**2. )

    dU_msh = get_dU(xe,re,R,CT,TI,order=1,pars=pars)
    dUeq = np.sum(weight_msh*dU_msh,axis=0)

    return dUeq

def GCLarsen_v0(WF, WS, WD, TI,
    pars=[0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0]):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Inputs
    ----------
    WF: WindFarm
        Windfarm instance
    WS: float
        Rotor averaged Undisturbed wind speed [m/s] for each WT
    WD: float
        Rotor averaged Undisturbed wind direction [deg] for each WT
        Meteorological axis. North = 0 [deg], clockwise.
    TI: float
        Rotor averaged turbulence intensity [-] for each WT

    Returns
    -------
    P_WT: ndarray
        Power production of the wind turbines (nWT,1) [W]
    U_WT: ndarray
        Wind speed at hub height (nWT,1) [m/s]
    Ct: ndarray
        Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    Dist, nDownstream, id0 = WF.turbineDistance(np.mean(WD))
    zg = WF.vectWTtoWT[2,:,:]

    # Initialize arrays to NaN
    Ct = np.nan * np.ones([WF.nWT])
    U_WT = copy(WS) #np.nan * np.ones([WF.nWT])
    P_WT = np.nan * np.ones([WF.nWT])

    # Initialize first upstream turbine
    Ct[id0[0]] = WF.WT[id0[0]].get_CT(WS[id0[0]])
    P_WT[id0[0]] = WF.WT[id0[0]].get_P(WS[id0[0]])
    U_WT[id0[0]] = WS[id0[0]]

    for i in range(1, WF.nWT):
        cWT = id0[i]  # Current wind turbine (wake operating)
        cR = WF.WT[cWT].R
        LocalDU = np.zeros([WF.nWT, 1])
        for j in range(i-1, -1, -1):
            # Loop on the upstream turbines (uWT) of the cWT
            uWT = id0[j]
            uWS = U_WT[uWT]  # Wind speed at wind turbine uWT
            uR = WF.WT[uWT].R
            uCT = Ct[uWT]
            if  np.isnan(uCT):
                uCT = WF.WT[uWT].get_CT(uWS)

            # WT2WT vector in wake coordinates
            Dist, _,_ = WF.turbineDistance(WD[uWT])
            x = Dist[0, uWT, cWT]
            y = Dist[1, uWT, cWT]
            z = zg[uWT, cWT]
            r = np.sqrt(y**2.+z**2.)

            # Calculate the wake width of uWT at the position of cWT
            Rw = get_Rw(x, uR, TI[uWT], uCT, pars)[0]
            if (r <= Rw+cR or uWS > 0):
                LocalDU[uWT] = uWS*get_dUeq(x,y,z,cR,uR,uCT,TI[uWT],pars)

        # Wake superposition
        DU = LocalDU.sum()

        U_WT[cWT] = U_WT[cWT] + DU
        if  U_WT[cWT] > WF.WT[cWT].u_cutin:
            Ct[cWT] = WF.WT[cWT].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[cWT].get_P(U_WT[cWT])
        else:
            Ct[cWT] = WF.WT[cWT].CT_idle
            P_WT[cWT] = 0.0

    return (P_WT,U_WT,Ct)

def GCLarsen(WF, WS, WD,TI,
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Inputs
    ----------
    WF: WindFarm
        Windfarm instance
    WS: float
        Undisturbed wind speed at hub height [m/s]
    WD: float
        Undisturbed wind direction at hub height [deg].
                Meteorological axis. North = 0 [deg], clockwise.
    TI: float
        Ambient turbulence intensity [-]

    Returns
    -------
    P_WT: ndarray
         Power production of the wind turbines (nWT,1) [W]
    U_WT: ndarray
         Wind speed at hub height (nWT,1) [m/s]
    Ct: float
        Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    distFlowCoord, nDownstream, id0 = WF.turbineDistance(np.mean(WD))
    zg = WF.vectWTtoWT[2,:,:]

    # Initialize arrays to NaN
    Ct = np.nan*np.ones([WF.nWT])
    P_WT = np.nan*np.ones([WF.nWT])

    # Initialize velocity to undisturbed eq ws
    U_WT  =  copy(WS)
    allR = np.array([WF.WT[i].R for i in range(WF.nWT)])

    # Extreme wake to define WT's in each wake, including partial wakes
    ID_wake = {i:((get_Rw(x=distFlowCoord[0,i,:],   # streamwise distance
                         R=allR[i],                # Upstream radius
                         CT=0.9,                   # Maximum effect
                         TI=0.01,
                         pars=pars)[0])*3.
                    > np.abs(distFlowCoord[1,i,:]) + allR).nonzero()[0]
               for i in id0}

    for i in range(WF.nWT):
        #Current wind turbine starting from the most upstream
        cWT = id0[i]
        # Current radius
        cR = WF.WT[cWT].R
        # Current wind speed
        cU = U_WT[cWT]
        if cU>WF.WT[cWT].u_cutin:
            Ct[cWT] = WF.WT[cWT].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[cWT].get_P(U_WT[cWT])
        else:
           Ct[cWT] = WF.WT[cWT].CT_idle # Drag coefficient of the idled turbine
           P_WT[cWT] = 0.0
        # Current turbine CT
        cCT=Ct[cWT]

        #Radial coordinates in cWT for wake affected WT's
        distFlowCoord,_,_ = WF.turbineDistance(WD[cWT])
        x = distFlowCoord[0, cWT, ID_wake[cWT]]
        y = distFlowCoord[1, cWT, ID_wake[cWT]]
        z = zg[cWT, ID_wake[cWT]]

        # Get all the rotor average wake deficits at the position of the -in wake-
        # downstream turbines
        localDU = cU*get_dUeq(x,y,z,allR[ID_wake[cWT]],cR,cCT,TI[cWT],pars)

        # Wake superposition
        U_WT[ID_wake[cWT]] = U_WT[ID_wake[cWT]] + localDU

        # Update Power and CT
        cU = U_WT[cWT]
        if cU>WF.WT[cWT].u_cutin:
            Ct[cWT] = WF.WT[cWT].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[cWT].get_P(U_WT[cWT])
        else:
           Ct[cWT] = WF.WT[cWT].CT_idle # Drag coefficient of the idled turbine
           P_WT[cWT] = 0.0

    return (P_WT,U_WT,Ct)

'''
print get_r96(D=80.0,CT=0.5,TI=0.05)
print
print get_Rw(x=np.linspace(10.,100.0,10),R=80.0,TI=0.05,a=0.3)
print
print get_dU(x=np.linspace(10.,100.0,2),r=np.linspace(10.,100.0,2).T,Rw=100.,U=10.0,R=80.0,TI=0.05,a=0.3)
print
print '----'
print get_dU(x=10., r=1., Rw=100., U=8., R=80., TI=0.07, a = 0.3)
print get_dU(x=10., r=101., Rw=100., U=8., R=80., TI=0.07, a = 0.3)
print '----'
print
print Ua(r=np.array([0.0,10.0]),te=np.array([[10.0],[-10.0]]),zc=80.0,us=15.0,z0=1.0)
print
print gaussN(R=80.0, func=Ua, varargin=[80.0,15.0,1.0], NG=4)
print
print gaussN(R=80.0, func=dU4Gauss, varargin=[80.0,0.0,100.,100.0,8.0,80.,0.05,0.3], NG=4)

v80 = wt.WindTurbine('Vestas v80 2MW offshore','V80_2MW_offshore.dat',70,40)
HR1 = wf.WindFarm('Horns Rev 1','HR_coordinates.dat',v80)

P_WT,U_WT = GCLarsen(WS=8.0,z0=0.0001,TI=0.05,WD=270,WF=HR1,NG=4,sup='quad',pars=[0.5,0.9,-0.124807893,0.136821858,15.6298,1.0])
print P_WT

'''
