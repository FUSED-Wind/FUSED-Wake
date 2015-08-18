"""GC Larsen wake model applied to offshore wind farms (WindFarm object)

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

References:
[1] Larsen GC. "A simple stationary semi-analytical wake model", 2009

"""
import numpy as np
import matplotlib.pyplot as plt
import WindTurbine as wt
import WindFarm as wf

def Ua(r,te,zc,us,z0):
    """Function of undisturbed inflow wind speed - log law.

    Inputs
    ----------
    r  (np.array float): Radial coord [m]
    te (np.array float): Angle coord where 0 is the horizontal right [rad]
    zc (float): Height to hub center [m]
    us (float): Friction velocity [m/s]
    z0 (float): Roughness height [m]

    Outputs
    ----------
    Ua  (np.array float): Axial wind speed [m/s]
    """
    kappa = 0.4 # Kappa: von karman constant
    return us / kappa * np.log((zc + r * np.sin(te)) / z0)

def Ua_shear(r,te,zc,uH,alpha):
    """Function of undisturbed inflow wind speed - power law.

    Inputs
    ----------
    r  (np.array float): Radial coord [m]
    te (np.array float): Angle coord where 0 is the horizontal right [rad]
    zc (float): Height to hub center [m]
    uH (float): Wind Speed at Hub Height [m/s]
    alpha (float): Shear Coefficient [-]

    Outputs
    ----------
    Ua  (np.array float): Axial wind speed [m/s]
    """

    return uH * ((zc + r * np.sin(te)) / zc)**alpha


def gaussN(R, func, varargin, NG=4):
    """Calculate numerically the gauss integration.
    [1] eq. 38

    Inputs
    ----------
    R (float): Wind turbine radius [m]
    func (function): Wind speed function
    varargin: Other arguments for the function besides [r,te]
    NG (int): Number of Ga

    Outputs
    ----------
    Ua (float):
    """
    A = np.pi*R**2
    #coefficients
    if  NG==4: #for speed give the full values
        rt = np.array([[ -0.339981043584856, -0.861136311594053,
            0.339981043584856, 0.861136311594053]])
        te = rt.T
        w = np.array([[0.652145154862546, 0.347854845137454,
            0.652145154862546, 0.347854845137454]])
    else:
        rt,w = np.polynomial.legendre.leggauss(NG)
        rt = np.array([rt])
        te = rt.T
        w = np.array([w])

    return (np.pi/4.0)*(R**2./A)*w*w.T*func(R*(rt+1.0)/2.0,
        np.pi*(te+1.0),*varargin)*(rt+1.0)

def get_R96(R,CT,TI,pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake radius at 9.6D downstream location of a turbine

    Inputs
    ----------
    R (float): Wind turbine radius
    CT (float): Outputs WindTurbine object's thrust coefficient
    TI (float): Ambient turbulence intensity

    Outputs
    ----------
    R96 (float): Wake radius at 9.6D downstream location
    """
    D = 2.0*R

    a1=pars[0]
    a2=pars[1]
    a3=pars[2]
    a4=pars[3]
    b1=pars[4]
    b2=pars[5]
    R96=a1*(np.exp(a2*CT*CT+a3*CT+a4*CT**0.))*(b1*TI+b2)*D

    return R96

def get_Rw(x,R,TI,CT,pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake radius at a location

    Inputs
    ----------
    x (float): Distance between turbines and wake location in the wind direction
    R (float): Wind turbine radius
    TI (float): Ambient turbulence intensity
    CT (float): Outputs WindTurbine object's thrust coefficient

    Outputs
    ----------
    Rw (float): Wake radius at a location
    """
    _ones = np.ones(np.shape(x))
    D = 2.0*R
    Area=np.pi*D*D/4.0

    #CT=4.0*a*(1.-a)

    m=1.0/(np.sqrt(1.0-CT))
    k=np.sqrt((m+1.0)/2.0)

    R96 = get_R96(R,CT,TI,pars)

    x0=(9.6*D)/((2.0*R96/(k*D))**3.0-1.0)
    term1=(k*D/2.0)**2.5
    term2=(105.0/(2.0*np.pi))**-0.5
    term3=(CT*Area*x0)**(-5.0/6.0)
    c1=term1*term2*term3

    Rw=((105.0*c1*c1/(2.0*np.pi))**0.2)*(CT*Area*(x+x0*_ones))**(1.0/3.0)

    if type(x)==float and x+x0 <=0. : Rw = 0
    elif type(x)==np.ndarray: Rw[x+x0*_ones<=0.] = 0.
    return Rw

def get_dU(x,r,Rw,U,R,TI,CT,
    order=1,
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake velocity deficit at a location

    Inputs
    ----------
    x (float): Distance between turbines and wake location in the wind direction
    r (float): Radial distance between the turbine and the location
    Rw (float): Wake radius at location [m]
    R (float): Wake producing turbine's radius [m]
    U (float): Undisturbed wind speed [m/s]
    TI (float): Ambient turbulence intensity [-]
    CT (float): Outputs WindTurbine object's thrust coefficient

    Outputs
    ----------
    dU (float): Wake velocity deficit at a location
    """
    _ones = np.ones(np.shape(x))

    D = 2.*R
    Area=np.pi*D*D/4.
    m=1./(np.sqrt(1.-CT))
    k=np.sqrt((m+1.)/2.)

    R96 = get_R96(R,CT,TI,pars)

    x0=(9.6*D)/((2.*R96/(k*D))**3.-1.)
    term1=(k*D/2.)**2.5
    term2=(105./(2.*np.pi))**-0.5
    term3=(CT*Area*x0)**(-5./6.)
    c1=term1*term2*term3

    term10=0.1111*U*_ones
    term20=(CT*Area*(x+x0)**(-2.))**(1./3.)
    term310=(r**(1.5))
    term320=(3.*c1*c1*CT*Area*(x+x0))**(-0.5)
    term30=term310*term320
    term40=((35./(2.*np.pi))**(3./10.))*((3.*c1*c1)**(-0.2))
    dU1=-term10*term20*(term30-term40)**2.

    dU = dU1

    if  order == 2:

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

        dU2_const = U*((CT*Area*((x+x0)**(-2.)))**(2./3.))
        dU2_term0 = d_0*(z**0.)
        dU2_term1 = d_1*(z**1.)
        dU2_term2 = d_2*(z**2.)
        dU2_term3 = d_3*(z**3.)
        dU2_term4 = d_4*(z**4.)

        dU2 = dU2_const*(dU2_term0+dU2_term1+dU2_term2+dU2_term3+dU2_term4)

        dU=dU1 + dU2

    if type(r)==np.ndarray: dU[Rw<r]=0. # Outside the wake
    elif type(r)==float and Rw<r: dU = 0

    if type(x)==np.ndarray: dU[x<=0.]=0. # upstream the wake gen. WT
    elif type(x)==float and x<=0.: dU = 0

    if CT==0: dU = 0.0*dU

    return dU

def dU4Gauss(r_e,t_e,rRi,tRi,x,Rw_,U,R,TI_,CT,
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    return get_dU(x=x,r=np.sqrt(rRi**2+r_e**2+2*rRi*r_e*np.cos(tRi-t_e)),
           Rw=Rw_,U=U,R=R,TI=TI_,CT=CT,pars=pars)

def GCLarsen_v0(WF,WS,WD,TI,z0,NG=4,sup='lin',
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Inputs
    ----------
    WF (WindFarm): Windfarm instance
    WS (float): Undisturbed wind speed at hub height [m/s]
    WD (float): Undisturbed wind direction at hub height [deg].
                Meteorological axis. North = 0 [deg], clockwise.
    TI (float): Ambient turbulence intensity [-]
    z0 (float): Roughness height [m]
    NG (int):   Number of points in Gaussian Quadrature for equivalent wind
                speed integration over rotor distFlowCoord
    sup (str):  Wake velocity deficit superposition method:
                'lin': Linear superposition
                'quad' Quadratic superposition

    Outputs
    ----------
    P_WT: Power production of the wind turbines (nWT,1) [W]
    U_WT: Wind speed at hub height (nWT,1) [m/s]
    Ct:   Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    (Dist,id0) = WF.turbineDistance(WD)

    kappa = 0.4 # Kappa: von karman constant
    us = WS*kappa/np.log(WF.WT.H/z0) #friction velocity
    WS_inf = gaussN(WF.WT.R, Ua, [WF.WT.H,us,z0]).sum() #eq inflow ws

    # Initialize arrays to NaN
    Ct = np.nan*np.ones([WF.nWT])
    U_WT = np.nan*np.ones([WF.nWT])
    P_WT = np.nan*np.ones([WF.nWT])
    MainWake = np.zeros([WF.nWT])
    MaxDU = 0.0

    # Initialize first upstream turbine
    Ct[id0[0]] = WF.WT.get_CT(WS_inf)
    U_WT[id0[0]] = WS_inf
    P_WT[id0[0]] = WF.WT.get_P(WS_inf)

    for i in range(1,WF.nWT):
        cWT = id0[i] #Current wind turbine
        cR = WF.WT.R
        LocalDU = np.zeros([WF.nWT,1])
        for j in range(i-1,-1,-1):
            # Loop on the upstream turbines of iWT
            uWT = id0[j]
            uWS = U_WT[uWT] # Wind speed at wind turbine uWT
            uCT = Ct[uWT]
            uR = WF.WT.R
            if  np.isnan(uCT):
                uCT = WF.WT.get_CT(uWS)

            WakeL = Dist[0,uWT,cWT]
            C2C =   Dist[1,uWT,cWT]

            # Calculate the wake width of jWT at the position of iWT
            Rw = get_Rw(WakeL,uR,TI,uCT,pars)
            if  (abs(C2C) <= Rw+cR or uWS==0):
                LocalDU[uWT] = gaussN(uR,dU4Gauss,[C2C,0.0,WakeL,Rw,uWS,uR,TI,uCT,pars],NG).sum()
                if  LocalDU[uWT]<MaxDU:
                    MaxDU = LocalDU[uWT]
                    MainWake[cWT]=uWT

        # Wake superposition
        if  sup=='lin':
            DU = LocalDU.sum()
        elif sup=='quad':
            DU = -np.sqrt(np.sum(LocalDU**2))

        U_WT[cWT] = max(0, WS_inf + DU)
        if  U_WT[cWT]>WF.WT.u_cutin:
            Ct[cWT] = WF.WT.get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT.get_P(U_WT[cWT])
        else:
           Ct[cWT] = 0.053
           P_WT[cWT] = 0.0

    return (P_WT,U_WT,Ct)

def GCLarsen(
    WF,
    WS,
    WD,
    TI,
    z0=0.0001,
    alpha=0.101,
    inflow='log',
    NG=4,
    sup='lin',
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Inputs
    ----------
    WF (WindFarm): Windfarm instance
    WS (float): Undisturbed wind speed at hub height [m/s]
    WD (float): Undisturbed wind direction at hub height [deg].
                Meteorological axis. North = 0 [deg], clockwise.
    TI (float): Ambient turbulence intensity [-]
    z0 (float): Roughness height [m]
    alpha (float): Shear coefficient [-]
                   Only used for power-law undisturbed inflow.
    inflow (Str):  Undisturbed inflow vertical profile:
                   'log': Logarithmic law (neutral case); uses z0
                   'pow': Power law profile; uses alpha
    NG (int):   Number of points in Gaussian Quadrature for equivalent wind
                speed integration over rotor distFlowCoord
    sup (str):  Wake velocity deficit superposition method:
                'lin': Linear superposition
                'quad' Quadratic superposition

    Outputs
    ----------
    P_WT: Power production of the wind turbines (nWT,1) [W]
    U_WT: Wind speed at hub height (nWT,1) [m/s]
    Ct:   Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    (distFlowCoord,id0) = WF.turbineDistance(WD)

    if inflow == 'log':
        kappa = 0.4 # Kappa: von karman constant
        us = WS*kappa/np.log(WF.WT[0].H/z0) #friction velocity
        #eq inflow ws
        WS_inf = gaussN(WF.WT[0].R, Ua, [WF.WT[0].H,us,z0]).sum()
    elif inflow == 'pow':
        #eq inflow ws
        WS_inf = gaussN(WF.WT[0].R, Ua_shear, [WF.WT[0].H,WS,alpha]).sum()

    # Gauss quadrature points
    r_Gc,w_Gc = np.polynomial.legendre.leggauss(NG)
    wj,wk=np.meshgrid(w_Gc,w_Gc)
    tj,rk=np.meshgrid(r_Gc,r_Gc)
    wj = wj.reshape((NG**2))
    tj = tj.reshape((NG**2))
    wk = wk.reshape((NG**2))
    rk = rk.reshape((NG**2))

    # Initialize arrays to NaN
    Ct = np.nan*np.ones([WF.nWT])
    P_WT = np.nan*np.ones([WF.nWT])

    # Initialize velocity to undisturbed eq ws
    U_WT  = WS_inf*np.ones([WF.nWT])
    U_WT0 = WS_inf*np.ones([WF.nWT])
    DU_sq = 0.*U_WT

    # Extreme wake to define WT's in each wake, including partial wakes
    ID_wake = {id0[i]:(get_Rw(x=distFlowCoord[0,id0[i],:],\
                              R=2.*WF.WT[0].R,TI=TI,CT=0.9,pars=pars)>\
               np.abs(distFlowCoord[1,id0[i],:])).nonzero()[0] \
               for i in range(WF.nWT)}

    for i in range(WF.nWT):
        #Current wind turbine starting from the most upstream
        cWT = id0[i]
        cR = WF.WT[i].R
        cU = U_WT[cWT]
        if  cU>WF.WT[i].u_cutin:
            Ct[cWT] = WF.WT[i].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[i].get_P(U_WT[cWT])
        else:
           Ct[cWT] = 0.053
           P_WT[cWT] = 0.0

        cCT=Ct[cWT]
        #ID_wake = {cWT:(get_Rw(x=distFlowCoord[0,cWT,:],\
        #                           R=cR*1.5,TI=TI,CT=cCT)\
        #           >np.abs(distFlowCoord[1,cWT,:])).nonzero()}

        #Radial coordinates in cWT for wake affected WT's
        x=distFlowCoord[0,cWT,ID_wake[cWT]]
        r_Ri  = np.abs(distFlowCoord[1,cWT,ID_wake[cWT]])
        th_Ri = np.pi*(np.sign(distFlowCoord[1,cWT,ID_wake[cWT]])+1.0)
        RW = get_Rw(x=x,R=cR,TI=TI,CT=cCT,pars=pars)

        # Meshgrids (Tensorial) extension of points of evaluation
        # to perform Gaussian quadrature
        r_Ri_m, rk_m = np.meshgrid(r_Ri, rk)
        th_Ri_m, tj_m = np.meshgrid(th_Ri, tj)
        x_m, wj_m = np.meshgrid(x, wj)
        RW_m, wk_m = np.meshgrid(RW, wk)

        # Radial points of evaluation
        r_eval = np.sqrt(r_Ri_m**2.0 + \
                 (cR*(rk_m+1.)/2.0)**2. + \
                 r_Ri_m*cR*(rk_m+1.)*\
                 np.cos(th_Ri_m-np.pi*(tj_m+1.)))

        # Eval wake velocity deficit
        DU_m = get_dU(x=x_m,r=r_eval,Rw=RW_m,
                      U=cU,R=cR,TI=TI,CT=cCT,pars=pars)

        localDU = np.sum((1./4.)*wj_m*wk_m*DU_m*(rk_m+1.0),axis=0)

        # Wake superposition
        if  sup=='lin':
            U_WT[ID_wake[cWT]] = U_WT[ID_wake[cWT]] + localDU
            U_WT[U_WT<0.]=0.
        elif sup=='quad':
            DU_sq[ID_wake[cWT]] = DU_sq[ID_wake[cWT]] + localDU**2.
            U_WT = U_WT0 - np.sqrt(DU_sq)
            U_WT[U_WT<0.]=0.

    return (P_WT,U_WT,Ct)

def GCL_P_GaussQ_Norm_U_WD(WF,WS,meanWD,stdWD,NG_P,TI,
    z0=0.0001,alpha=0.101,inflow='log',NG=4,sup='lin',
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the Gaussian quadrature average of GCLarsen
    power/WS prediction under normally distributed wind direction
    uncertainty inside the Reynolds averaging time

    Inputs
    ----------
    meanWD (float): Mean wind direction [deg]
    stdWD (float): Std wind direction [deg]
    NG_P (int): Number of Gaussian quadrature points for power avg.
    WS (float): Undisturbed wind speed at hub height [m/s]
    z0 (float): Roughness height [m]
    TI (float): Ambient turbulence intensity

    Outputs
    ----------
    P_WT: Mean Power production of the wind turbines (nWT,1) [W]
    U_WT: Mean Wind speed at hub height (nWT,1) [m/s]
    CT_WT: Mean thrust coefficient [-]
    """

    xi, wi = np.polynomial.hermite.hermgauss(NG_P)

    meanP_WT = np.zeros([WF.nWT])
    meanU_WT = np.zeros([WF.nWT])
    meanCT_WT = np.zeros([WF.nWT])
    for i in range(NG_P):
        P_WT,U_WT,CT_WT = GCLarsen(
                            WF=WF,
                            WS=WS,
                            WD=meanWD+np.sqrt(2.)*stdWD*xi[i],
                            TI=TI,
                            z0=z0,
                            alpha=alpha,
                            inflow=inflow,
                            NG=NG,
                            sup=sup,
                            pars=pars)
        meanP_WT += wi[i]*P_WT*(1./np.sqrt(np.pi))
        meanU_WT += wi[i]*U_WT*(1./np.sqrt(np.pi))
        meanCT_WT += wi[i]*CT_WT*(1./np.sqrt(np.pi))

    return meanP_WT,meanU_WT,meanCT_WT

def GCL_P_GaussQ_Uni_U_WD(WF,WS,meanWD,U_WD,NG_P,TI,
    z0=0.0001,alpha=0.101,inflow='log',NG=4,sup='lin',
    pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the Gaussian quadrature average of GCLarsen
    power/WS prediction under normally distributed wind direction
    uncertainty

    Inputs
    ----------
    meanWD (float): Mean wind direction [deg]
    stdWD (float): Std wind direction [deg]
    NG_P (int): Number of Gaussian quadrature points for power avg.
    WS (float): Undisturbed wind speed at hub height [m/s]
    z0 (float): Roughness height [m]
    TI (float): Ambient turbulence intensity

    Outputs
    ----------
    P_WT: Mean Power production of the wind turbines (nWT,1) [W]
    U_WT: Mean Wind speed at hub height (nWT,1) [m/s]
    CT_WT: Mean thrust coefficient [-]
    """

    xi, wi = np.polynomial.legendre.leggauss(NG_P)

    meanP_WT = np.zeros([WF.nWT])
    meanU_WT = np.zeros([WF.nWT])
    meanCT_WT = np.zeros([WF.nWT])
    for i in range(NG_P):
        P_WT,U_WT,CT_WT = GCLarsen(
                            WF=WF,
                            WS=WS,
                            WD=meanWD+U_WD*xi[i],
                            TI=TI,
                            z0=z0,
                            alpha=alpha,
                            inflow=inflow,
                            NG=NG,
                            sup=sup,
                            pars=pars)
        meanP_WT += wi[i]*P_WT*1./2.
        meanU_WT += wi[i]*U_WT*1./2.
        meanCT_WT += wi[i]*CT_WT*1./2.

    return meanP_WT,meanU_WT,meanCT_WT
'''
print get_R96(R=80.0,CT=0.5,TI=0.05)
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
