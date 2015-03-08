"""An offshore wind farm model

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

"""
import numpy as np
import matplotlib.pyplot as plt
import WindTurbine as wt

class WindFarm:
    def __init__(self, name, coordFile, WT):
        """Initializes a WindFarm object

        Inputs
        ----------
        name (str):  WindFarm name
        coordFile (str): Wind Farm layout coordinates text file.
        WindTurbine (WindTurbine): WindTurbine object (only one type per WindFarm)

        Outputs
        ----------
        WindTurbine (WindFarm)
        """
        self.name = name
        coordArray = np.loadtxt(coordFile)
        self.pos = coordArray.T # np.array(2 x nWT)
        self.nWT = len(self.pos[0,:])
        self.WT = WT

        # Vector from iWT to jWT: self.vectWTtoWT[:,i,j]
        self.vectWTtoWT=np.swapaxes([self.pos-np.repeat(np.atleast_2d \
            (self.pos[:,i]).T,self.nWT,axis=1) for i in range(self.nWT)],0,1)

    def __repr__(self):
        print "-------------------------------------------------------------"
        print "%s has %s %s wind turbines " %(self.name,self.nWT,self.WT.name)
        print "-------------------------------------------------------------"
        return ''

    def turbineDistance(self,wd):
        """Computes the WT to WT distance in flow coordinates
        ranks the most of most upstream turbines

        Inputs
        ----------
        wd (float): Wind direction in degrees

        Outputs
        ----------
        distFlowCoord: Vector from iWT to jWT: self.vectWTtoWT[:,i,j]
        idWT (int np.array): turbine index array
        """
        angle = np.radians(270.-wd)
        ROT = np.array([[np.cos(angle), np.sin(angle)], \
              [-np.sin(angle), np.cos(angle)]])
        distFlowCoord = np.einsum('ij,jkl->ikl',ROT,self.vectWTtoWT)
        nDownstream = [(distFlowCoord[0,i,:]<0).sum() \
                      for i in range(self.nWT)]
        ID0 = np.argsort(nDownstream)
        return distFlowCoord, ID0

    def toFlowCoord(self,wd,vect):
        """Rotates a 2xN np.array to flow coordinates

        Inputs
        ----------
        wd (float): Wind direction in degrees
        vect (np.array): Vector or Matrix 2xN

        Outputs
        ----------
        vect (np.array): Vector or Matrix 2xN
        """
        angle = np.radians(270.-wd)
        ROT = np.array([[np.cos(angle), np.sin(angle)],\
              [-np.sin(angle), np.cos(angle)]])
        return np.dot(ROT,vect)

    def plot(self,WT_num=False):
        x = (self.pos[0,:]-min(self.pos[0,:]))/(2.*self.WT.R)
        y = (self.pos[1,:]-min(self.pos[1,:]))/(2.*self.WT.R)
        fig, ax = plt.subplots()
        ax.scatter(x,y,c='black')
        if WT_num==True:
            for i in range(0,self.nWT):
                ax.annotate(i, (x[i],y[i]))
        elif WT_num==False:
            print 'No annotation of turbines'
        ax.set_xlabel('x/D [-]')
        ax.set_ylabel('y/D [-]')
        ax.axis('equal')
        ax.set_title(self.name)
        return fig,ax

    def plot_order(self,wd):
        x = (self.pos[0,:]-min(self.pos[0,:]))/1000
        y = (self.pos[1,:]-min(self.pos[1,:]))/1000
        dist,nDownstream,idWT = self.turbineDistance(wd)
        fig, ax = plt.subplots()
        ax.scatter(x,y,c='black')
        for i in range(0,self.nWT):
            ax.annotate(int(idWT[i]), (x[i],y[i]))
        ax.set_xlabel('x [km]')
        ax.set_ylabel('y [km]')
        ax.set_title(self.name+' Wind direction '+str(wd))
        return fig,ax

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