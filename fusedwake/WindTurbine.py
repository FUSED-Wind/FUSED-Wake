"""An offshore wind farm model

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

"""
import numpy as np
try:
    import scipy as sp
    from scipy.interpolate import interp1d as interpolator
except Exception as e:
    print(e.message)
# from scipy.interpolate import pchipInterpolator as interpolator

class WindTurbine:
    """Wind Turbine instance

    Defines WT parameters and operational curves

    """
    def __init__(self,name,refCurvesFile,H,R,CT_idle=0.053):
        """Initializes a WindTurbine object

        Parameters
        ----------
        name: str
            Wind turbine name
        refCurvesFile: str
            Power and thrust coefficient curves text file.
        H: float
            Wind turbine hub height [m]
        R: float
            Radius [m]

        Returns
        -------
        WindTurbine (WindTurbine)
        """
        self.name = name
        self.H = H
        self.R = R

        refCurvesArray=np.loadtxt(refCurvesFile,delimiter=', ',skiprows=5)

        self.refCurvesArray = refCurvesArray

        self.CT_idle = CT_idle

        self.wt_init()

    def wt_init(self):
        self.ref_u = self.refCurvesArray[:,0]
        self.ref_P = self.refCurvesArray[:,1]
        self.ref_CT = self.refCurvesArray[:,2]
        self.u_cutin = self.ref_u[0]
        self.u_cutout = self.ref_u[-1]
        self.P_rated = np.max(self.ref_P)

        self.PCI = interpolator(self.ref_u, self.ref_P)
        self.CTCI = interpolator(self.ref_u, self.ref_CT)

        index = np.nonzero(self.ref_P==self.P_rated)[0][0]
        self.PCI_u = interpolator(
                     self.ref_P[:index+1],self.ref_u[:index+1])

        self.u_rated = np.float(self.PCI_u(self.P_rated))

    def __repr__(self):
        print("-------------------------------------------------------------")
        print("\t %s" % (self.name))
        print("-------------------------------------------------------------")
        print("\nHeight \t %s [m]\nRadius \t %s [m] \n" %(self.H, self.R))
        print("-------------------------------------------------------------")
        print("\t u [m/s] \t P [kW] \t CT [-]")
        for row in self.refCurvesArray:
            print('\t %0.0f \t\t %0.0f \t\t %0.3f'%(row[0],row[1]/1000.0,row[2]))
        print("-------------------------------------------------------------")
        return ''


    def get_P(self,u):
        """Computes the Power of the WindTurbine
        at the undisturbed wind speed

        Parameters
        ----------
        u: float
            Undisturbed wind speed

        Returns
        -------
        Power: float
            WindTurbine object's power
        """
        #return np.interp(u, self.ref_u, self.ref_P, left=0)
        return ((u>=self.u_cutin)&(u<=self.u_cutout))*self.PCI(u)

    def get_u(self, P):
        """Computes the undisturbed wind speed of the WindTurbine
        given a power under rated power

        Parameters
        ----------
        P: float
            Power

        Returns
        -------
        u: float
            Undisturbed wind speed
        """
        return ((P >= 0.0) & (P <= self.P_rated)) * self.PCI_u(P)

    def get_CT(self,u):
        """Computes the thrust coefficient of the WindTurbine
        at the undisturbed wind speed

        Parameters
        ----------
        u: float
            Undisturbed wind speed

        Returns
        -------
        CT: float
            Thrust coefficient
        """
        #return np.interp(u, self.ref_u, self.ref_CT)
        return ((u>=self.u_cutin)&(u<=self.u_cutout))*self.CTCI(u) + \
                ((u<self.u_cutin)|(u>self.u_cutout))*self.CT_idle

    def get_a(self,CT):
        """Computes the axial velocity deficit coefficient at the rotor disc
        bosed on the WindTurbine's thrust coefficient

        .. math::
            a = \\frac{1}{2} \\left( 1 - \\sqrt{1-C_T} \\right)

        Parameters
        ----------
        CT: float
            Thrust coefficient

        Returns
        -------
        a: float
            Axial velocity deficit coefficient at the rotor disc
        """
        return 0.5 * ( 1. - np.sqrt(1.-CT))


'''
v80 = WindTurbine('Vestas v80 2MW offshore','V80_2MW_offshore.dat',70,40)
v80.display_windTurbine()
print(v80.get_P(u=np.array([10.5])))
print(v80.get_CT(u=np.array([10.5])))
'''
