
"""GC Larsen wake model applied to offshore wind farms (WindFarm object)
@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>
References:
[1] Larsen GC. "A simple stationary semi-analytical wake model", 2009
"""

import numpy as np

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
