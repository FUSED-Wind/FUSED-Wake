"""An offshore wind farm model

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

"""
import numpy as np
MATPLOTLIB = True
try:
    import matplotlib.pyplot as plt
except Exception as e:
    MATPLOTLIB = False
    print("WARNING: Matplotlib isn't installed correctly:", e)

class WindTurbineList(list):
    """A simple list class that can also act as a single element when needed.
    Accessing one of the attribute of this list will get the first element of the
    list attributes.

    Note
    ----
    We assume here that if a model call this intense as if it's one wind
    turbine, the user has been clever enough to pass as input a wind farm with
    identical turbines.
    """
    def __getattr__(self, key):
        #TODO: make some checks to catch possible bugs when the turbines are not similar.
        return getattr(self.__getitem__(0), key)

    def names(self):
        return [getattr(w, 'name') for w in self]


class WindFarm(object):
    def __init__(self, name, coordFile, WT):
        """Initializes a WindFarm object

        Parameters
        ----------
        name: str
            WindFarm name
        coordFile: str
            Wind Farm layout coordinates text file.
        WindTurbine: WindTurbine
            WindTurbine object (only one type per WindFarm)
        """
        self.name = name
        coordArray = np.loadtxt(coordFile)
        self.pos = coordArray.T  # np.array(2 x nWT)
        self.nWT = len(self.pos[0, :])
        # We generate a wind turbine list
        self.WT = WindTurbineList([WT for i in range(self.nWT)])

        # Vector from iWT to jWT: self.vectWTtoWT[:,i,j]
        self.vectWTtoWT = np.swapaxes([self.pos -
            np.repeat(np.atleast_2d(self.pos[:, i]).T, self.nWT, axis=1)
            for i in range(self.nWT)], 0, 1)


    def __repr__(self):
        print("-------------------------------------------------------------")
        print("%s has %s %s wind turbines "%(self.name, self.nWT, self.WT.name))
        print("-------------------------------------------------------------")
        return ''

    def turbineDistance(self, wd):
        """Computes the WT to WT distance in flow coordinates
        ranks the most of most upstream turbines

        Parameters
        ----------
        wd: float
            Wind direction in degrees

        Returns
        -------
        distFlowCoord: Vector from iWT to jWT: self.vectWTtoWT[:,i,j]
        idWT: ndarray(int)
            turbine index array
        """
        angle = np.radians(270.-wd)
        ROT = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])
        distFlowCoord = np.einsum('ij,jkl->ikl', ROT, self.vectWTtoWT)
        nDownstream = [(distFlowCoord[0, i, :] < 0).sum() for i in range(self.nWT)]
        ID0 = np.argsort(nDownstream)
        return distFlowCoord, nDownstream, ID0

    def toFlowCoord(self, wd, vect):
        """Rotates a 2xN np.array to flow coordinates

        Parameters
        ----------
        wd: float
            Wind direction in degrees
        vect: ndarray
            Vector or Matrix 2xN

        Returns
        -------
        vect: ndarray
            Vector or Matrix 2xN
        """
        angle = np.radians(270.-wd)
        ROT = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])
        return np.dot(ROT, vect)

    def plot(self, WT_num=False):
        if MATPLOTLIB:
            x = (self.pos[0, :] - min(self.pos[0, :])) / (2. * self.WT.R)
            y = (self.pos[1, :] - min(self.pos[1, :])) / (2. * self.WT.R)
            fig, ax = plt.subplots()
            ax.scatter(x, y, c='black')
            if WT_num:
                for i in range(0, self.nWT):
                    ax.annotate(i, (x[i], y[i]))
            elif not WT_num:
                print('No annotation of turbines')
            ax.set_xlabel('x/D [-]')
            ax.set_ylabel('y/D [-]')
            ax.axis('equal')
            ax.set_title(self.name)
            return fig, ax

    def plot_order(self, wd):
        x = (self.pos[0, :] - min(self.pos[0, :])) / 1000
        y = (self.pos[1, :] - min(self.pos[1, :])) / 1000
        dist, nDownstream, idWT = self.turbineDistance(wd)
        fig, ax = plt.subplots()
        ax.scatter(x, y, c='black')
        for i in range(0, self.nWT):
            ax.annotate(int(idWT[i]), (x[i], y[i]))
        ax.set_xlabel('x [km]')
        ax.set_ylabel('y [km]')
        ax.set_title(self.name+' Wind direction '+str(wd))
        return fig, ax

    def __getattr__(self, key):
        """Give access to a list of the properties of the turbine

        Parameters
        ----------
        key: str
            The parameter to return

        Returns
        -------
        parameters: list
            The parameter list of the turbines
        """
        return [getattr(wt, key) for wt in self.WT]

'''
v80 = wt.WindTurbine('Vestas v80 2MW offshore','V80_2MW_offshore.dat',70,40)
HR1 = WindFarm('Horns Rev 1','HR_coordinates.dat',v80)
HR1.display_windFarm()

HR1.plot()
HR1.plot_order(0)
HR1.plot_order(90)
HR1.plot_order(180)
HR1.plot_order(270)
'''
