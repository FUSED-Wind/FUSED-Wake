
if inflow == 'log':
    kappa = 0.4 # Kappa: von karman constant
    us = WS * kappa / np.log(WF.WT[0].H / z0)  # friction velocity
    #eq inflow ws
    WS_inf = gaussN(WF.WT[0].R, Ua, [WF.WT[0].H,us,z0]).sum()
elif inflow == 'pow':
    #eq inflow ws
    WS_inf = gaussN(WF.WT[0].R, Ua_shear, [WF.WT[0].H,WS,alpha]).sum()


def gauss_rws()

class WFFM(WindFarm):
    def supperposition(self, ws, **kwargs):
        """
        """
        pass

    def single_wake(self, x, r, **kwargs):
        """
        """
        pass

    def wake_width(self, x, **kwargs):
        """
        """
        pass

    def rotor_wind_speed(self, ws, **kwargs):
        """
        """
        pass

    def run(self):
        """
        """
        (distFlowCoord,id0) = WF.turbineDistance(WD)

        # TODO: decide how at what height the us should be defined
        if inflow == 'log':
            kappa = 0.4 # Kappa: von karman constant
            us = WS * kappa / np.log(WF.WT[0].H / z0)  # friction velocity
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

        allR = np.array([WF.WT[i].R for i in range(WF.nWT)])

        # Extreme wake to define WT's in each wake, including partial wakes
        ID_wake = {i:(get_Rw(x=distFlowCoord[0,i,:],                # streamwise distance
                             R=WF.WT[i].R,                          # Upstream radius
                             TI=TI,
                             CT=0.9,                                #Maximum effect
                             pars=pars)
                            > np.abs(distFlowCoord[1,i,:]) + allR).nonzero()[0]
                   for i in id0}

        for i in range(WF.nWT):
            #Current wind turbine starting from the most upstream
            cWT = id0[i]
            # Current radius
            cR = WF.WT[i].R
            # Current hub wind speed
            cU = U_WT[cWT]
            if cU>WF.WT[i].u_cutin:
                Ct[cWT] = WF.WT[i].get_CT(U_WT[cWT])
                P_WT[cWT] = WF.WT[i].get_P(U_WT[cWT])
            else:
               Ct[cWT] = 0.053  # Drag coefficient of the idled turbine
               P_WT[cWT] = 0.0

            # Current turbine CT
            cCT=Ct[cWT]
            #ID_wake = {cWT:(get_Rw(x=distFlowCoord[0,cWT,:],\
            #                           R=cR*1.5,TI=TI,CT=cCT)\
            #           >np.abs(distFlowCoord[1,cWT,:])).nonzero()}

            #Radial coordinates in cWT for wake affected WT's
            x = distFlowCoord[0, cWT, ID_wake[cWT]]
            r_Ri  = np.abs(distFlowCoord[1,cWT, ID_wake[cWT]])
            th_Ri = np.pi*(np.sign(distFlowCoord[1, cWT, ID_wake[cWT]]) + 1.0) # <- what is this? [0|2pi]

            # Get all the wake radius at the position of the -in wake- downstream turbines
            RW = get_Rw(x=x, R=cR, TI=TI, CT=cCT, pars=pars)

            # Meshgrids (Tensorial) extension of points of evaluation
            # to perform Gaussian quadrature
            r_Ri_m, rk_m = np.meshgrid(r_Ri, rk)
            th_Ri_m, tj_m = np.meshgrid(th_Ri, tj)
            x_m, wj_m = np.meshgrid(x, wj)
            RW_m, wk_m = np.meshgrid(RW, wk)

            # downstream Radius
            downR = np.array([WF.WT[j].R for j in ID_wake[cWT]])
            downR_m, dummyvar = np.meshgrid(downR, np.zeros((NG**2)))
            cH = WF.WT[cWT].H
            downH = np.array([WF.WT[j].H for j in ID_wake[cWT]]) - cH
            downH_m, dummyvar = np.meshgrid(downH, np.zeros((NG**2)))

            # Radial points of evaluation    <- probably need to add the turbine height difference here?
            r_eval = np.sqrt(r_Ri_m**2.0 +
                             (downR_m * (rk_m + 1.) / 2.0)**2. +
                             r_Ri_m * downR_m * (rk_m + 1.) * np.cos(th_Ri_m - np.pi*(tj_m + 1.)))

            # Eval wake velocity deficit
            DU_m = get_dU(x=x_m, r=r_eval, Rw=RW_m,
                          U=cU, R=downR_m, TI=TI, CT=cCT, pars=pars)

            # Gaussian average
            localDU = np.sum((1./4.) * wj_m * wk_m * DU_m * (rk_m + 1.0), axis=0)

            # Wake superposition
            if sup == 'lin':
                U_WT[ID_wake[cWT]] = U_WT[ID_wake[cWT]] + localDU
                U_WT[U_WT<0.]=0.
            elif sup == 'quad':
                DU_sq[ID_wake[cWT]] = DU_sq[ID_wake[cWT]] + localDU**2.
                U_WT = U_WT0 - np.sqrt(DU_sq)
                U_WT[U_WT<0.]=0.

        return (P_WT,U_WT,Ct)

    def inflow_wind_speed(self):

    def __call__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.run()



def wffm_vectorized(WF, WS, WD,TI,
    z0=0.0001, alpha=0.101, inflow='log', NG=4, sup='lin'):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Parameters
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
    z0: float, optional
        Roughness height [m]
    alpha: float, optional
        Shear coefficient [-]
        Only used for power-law undisturbed inflow.
    inflow: Str, optional
        Undisturbed inflow vertical profile:
           'log': Logarithmic law (neutral case); uses z0
           'pow': Power law profile; uses alpha
    NG: int, optional
        Number of points in Gaussian Quadrature for equivalent wind
        speed integration over rotor distFlowCoord
    sup: str, optional
        Wake velocity deficit superposition method:
            'lin': Linear superposition
            'quad' Quadratic superposition

    Returns
    -------
    P_WT: ndarray
         Power production of the wind turbines (nWT,1) [W]
    U_WT: ndarray
         Wind speed at hub height (nWT,1) [m/s]
    Ct: float
        Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    (distFlowCoord,id0) = WF.turbineDistance(WD)

    # TODO: decide how at what height the us should be defined
    if inflow == 'log':
        kappa = 0.4 # Kappa: von karman constant
        us = WS * kappa / np.log(WF.WT[0].H / z0)  # friction velocity
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

    allR = np.array([WF.WT[i].R for i in range(WF.nWT)])

    # Extreme wake to define WT's in each wake, including partial wakes
    ID_wake = {i:(get_Rw(x=distFlowCoord[0,i,:],                # streamwise distance
                         R=WF.WT[i].R,                          # Upstream radius
                         TI=TI,
                         CT=0.9,                                #Maximum effect
                         pars=pars)
                        > np.abs(distFlowCoord[1,i,:]) + allR).nonzero()[0]
               for i in id0}

    for i in range(WF.nWT):
        #Current wind turbine starting from the most upstream
        cWT = id0[i]
        # Current radius
        cR = WF.WT[i].R
        # Current hub wind speed
        cU = U_WT[cWT]
        if cU>WF.WT[i].u_cutin:
            Ct[cWT] = WF.WT[i].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[i].get_P(U_WT[cWT])
        else:
           Ct[cWT] = 0.053  # Drag coefficient of the idled turbine
           P_WT[cWT] = 0.0

        # Current turbine CT
        cCT=Ct[cWT]
        #ID_wake = {cWT:(get_Rw(x=distFlowCoord[0,cWT,:],\
        #                           R=cR*1.5,TI=TI,CT=cCT)\
        #           >np.abs(distFlowCoord[1,cWT,:])).nonzero()}

        #Radial coordinates in cWT for wake affected WT's
        x = distFlowCoord[0, cWT, ID_wake[cWT]]
        r_Ri  = np.abs(distFlowCoord[1,cWT, ID_wake[cWT]])
        th_Ri = np.pi*(np.sign(distFlowCoord[1, cWT, ID_wake[cWT]]) + 1.0) # <- what is this? [0|2pi]

        # Get all the wake radius at the position of the -in wake- downstream turbines
        RW = get_Rw(x=x, R=cR, TI=TI, CT=cCT, pars=pars)

        # Meshgrids (Tensorial) extension of points of evaluation
        # to perform Gaussian quadrature
        r_Ri_m, rk_m = np.meshgrid(r_Ri, rk)
        th_Ri_m, tj_m = np.meshgrid(th_Ri, tj)
        x_m, wj_m = np.meshgrid(x, wj)
        RW_m, wk_m = np.meshgrid(RW, wk)

        # downstream Radius
        downR = np.array([WF.WT[j].R for j in ID_wake[cWT]])
        downR_m, dummyvar = np.meshgrid(downR, np.zeros((NG**2)))
        cH = WF.WT[cWT].H
        downH = np.array([WF.WT[j].H for j in ID_wake[cWT]]) - cH
        downH_m, dummyvar = np.meshgrid(downH, np.zeros((NG**2)))

        # Radial points of evaluation    <- probably need to add the turbine height difference here?
        r_eval = np.sqrt(r_Ri_m**2.0 +
                         (downR_m * (rk_m + 1.) / 2.0)**2. +
                         r_Ri_m * downR_m * (rk_m + 1.) * np.cos(th_Ri_m - np.pi*(tj_m + 1.)))

        # Eval wake velocity deficit
        DU_m = get_dU(x=x_m, r=r_eval, Rw=RW_m,
                      U=cU, R=downR_m, TI=TI, CT=cCT, pars=pars)

        # Gaussian average
        localDU = np.sum((1./4.) * wj_m * wk_m * DU_m * (rk_m + 1.0), axis=0)

        # Wake superposition
        if sup == 'lin':
            U_WT[ID_wake[cWT]] = U_WT[ID_wake[cWT]] + localDU
            U_WT[U_WT<0.]=0.
        elif sup == 'quad':
            DU_sq[ID_wake[cWT]] = DU_sq[ID_wake[cWT]] + localDU**2.
            U_WT = U_WT0 - np.sqrt(DU_sq)
            U_WT[U_WT<0.]=0.

    return (P_WT,U_WT,Ct)


def wffm_classic(WF, WS, WD, TI, z0, NG=4, sup='lin',
    pars=[0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0]):
    """Computes the WindFarm flow and Power using GCLarsen
    [Larsen, 2009, A simple Stationary...]

    Parameters
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
    z0: float
        Roughness height [m]

    NG: int, optional
        Number of points in Gaussian Quadrature for equivalent wind
        speed integration over rotor distFlowCoord
    sup: str, optional
        Wake velocity deficit superposition method:
            'lin': Linear superposition
            'quad' Quadratic superposition


    Returns
    -------
    P_WT: ndarray
        Power production of the wind turbines (nWT,1) [W]
    U_WT: ndarray
        Wind speed at hub height (nWT,1) [m/s]
    Ct: ndarray
        Thrust coefficients for each wind turbine (nWT,1) [-]
    """
    (Dist, id0) = WF.turbineDistance(WD)

    kappa = 0.4  # Kappa: von karman constant
    us = WS * kappa/np.log(WF.WT[id0[0]].H/z0)  # friction velocity
    WS_inf = gaussN(WF.WT[id0[0]].R, Ua, [WF.WT[id0[0]].H, us, z0]).sum()  # eq inflow ws

    # Initialize arrays to NaN
    Ct = np.nan * np.ones([WF.nWT])
    U_WT = np.nan * np.ones([WF.nWT])
    P_WT = np.nan * np.ones([WF.nWT])
    MainWake = np.zeros([WF.nWT])
    MaxDU = 0.0

    # Initialize first upstream turbine
    Ct[id0[0]] = WF.WT[id0[0]].get_CT(WS_inf)
    U_WT[id0[0]] = WS_inf
    P_WT[id0[0]] = WF.WT[id0[0]].get_P(WS_inf)

    for i in range(1, WF.nWT):
        cWT = id0[i]  # Current wind turbine
        cR = WF.WT[cWT].R
        LocalDU = np.zeros([WF.nWT, 1])
        for j in range(i-1, -1, -1):
            # Loop on the upstream turbines of iWT
            uWT = id0[j]
            uWS = U_WT[uWT]  # Wind speed at wind turbine uWT
            uCT = Ct[uWT]
            uR = WF.WT[uWT].R
            if np.isnan(uCT):
                uCT = WF.WT.get_CT(uWS)

            WakeL = Dist[0, uWT, cWT]
            C2C =   Dist[1, uWT, cWT]

            # Calculate the wake width of jWT at the position of iWT
            Rw = get_Rw(WakeL, uR, TI, uCT, pars)
            if (abs(C2C) <= Rw+cR or uWS == 0):
                LocalDU[uWT] = gaussN(uR, dU4Gauss,
                    [C2C, 0.0, WakeL, Rw, uWS, uR, TI, uCT, pars], NG).sum()
                if LocalDU[uWT]<MaxDU:
                    MaxDU = LocalDU[uWT]
                    MainWake[cWT] = uWT

        # Wake superposition
        if sup == 'lin':
            DU = LocalDU.sum()
        elif sup == 'quad':
            DU = -np.sqrt(np.sum(LocalDU**2))

        U_WT[cWT] = max(0, WS_inf + DU)
        if U_WT[cWT] > WF.WT[cWT].u_cutin:
            Ct[cWT] = WF.WT[cWT].get_CT(U_WT[cWT])
            P_WT[cWT] = WF.WT[cWT].get_P(U_WT[cWT])
        else:
            Ct[cWT] = 0.053
            P_WT[cWT] = 0.0

    return (P_WT,U_WT,Ct)
