import python.gcl as gcl
import fortran as fgcl
import numpy as np


class GCL(object):
    # The different versions and their respective inputs
    inputs = {
        'py0': ['WF', 'WS', 'WD', 'TI', 'z0', 'NG', 'sup', 'pars'],
        'py1': ['WF', 'WS', 'WD', 'TI', 'z0', 'alpha', 'inflow', 'NG', 'sup', 'pars'],
        'fort0': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ti',
                  'av', 'a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'ng', 'rho', 'ws_ci',
                  'ws_co', 'ct_idle'],
    }
    # Default variables for running the wind farm flow model
    defaults = {
        'rho': 1.225,
        'version': 'py1',
        'z0': 0.0001,
        'sup': 'lin', # ['lin' | 'quad']
        'alpha': 0.101,
        'pars': [0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0],
        'inflow': 'log',
        'NG': 4,
    }
    def __init__(self, **kwargs):
        self.set(self.defaults)
        self.set(kwargs)

    def set(self, dic):
        """ Set the attributes of a dictionary as instance variables. Prepares
        for the different versions of the wake model

        Parameters
        ----------
        dic: dict
            An input dictionary
        """
        for k, v in dic.items():
            setattr(self, k, v)

        # Preparing for the inputs for the fortran version
        if 'WF' in dic:
            self.x_g, self.y_g, self.z_g = self.WF.get_T2T_gl_coord2()
            self.dt = self.WF.rotor_diameter
            self.p_c = self.WF.power_curve
            self.ct_c = self.WF.c_t_curve
            self.ws_ci = self.WF.cut_in_wind_speed
            self.ws_co = self.WF.cut_out_wind_speed
            self.ct_idle = self.WF.c_t_idle


    def _get_kwargs(self, version):
        """Prepare a dictionary of inputs to be passed to the wind farm flow model

        Parameters
        ----------
        version: str
            The version of the wind farm flow model to run ['py0' | 'py1' | 'fort0']
        """
        if version == 'fort0':
            if isinstance(self.WS, float) or isinstance(self.WS, int):
                self.ws = np.array([self.WS])
                self.wd = np.array([self.WD])
                self.ti = np.array([self.TI])
            if not hasattr(self, 'wt_available'):
                self.wt_available = np.ones([len(self.ws), self.WF.nWT])
                self.av = self.wt_available
            elif self.wt_available.shape == (len(self.ws), self.WF.nWT):
                self.av = self.wt_available
            else:
                # stacking up the availability vector for each flow case
                self.av = np.vstack([self.wt_available for i in range(len(self.ws))])
            self.a1 = self.pars[0] * np.ones_like(self.ws)
            self.a2 = self.pars[1] * np.ones_like(self.ws)
            self.a3 = self.pars[2] * np.ones_like(self.ws)
            self.a4 = self.pars[3] * np.ones_like(self.ws)
            self.b1 = self.pars[4] * np.ones_like(self.ws)
            self.b2 = self.pars[5] * np.ones_like(self.ws)
            self.ng = self.NG

        return {k:getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}

    def __call__(self, *args, **kwargs):
        self.set(kwargs)
        if hasattr(self, 'version'):
            if   self.version == 'py0':
                self.p_wt, self.u_wt, self.c_t = gcl.GCLarsen_v0(**self._get_kwargs(self.version))
            elif self.version == 'py1':
                self.p_wt, self.u_wt, self.c_t = gcl.GCLarsen(**self._get_kwargs(self.version))
            elif self.version == 'fort0':
                #print(self._get_kwargs(self.version))
                try:
                    self.p_wt, self.t_wt, self.u_wt = fgcl.gcl_av(**self._get_kwargs(self.version))
                except Exception as e:
                    raise Exception('The fortran model failed with the followind inputs:', self._get_kwargs(self.version), e.msg)
                A = 0.25 * self.WF.WT.rotor_diameter**2.0
                self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
                self.p_wt *= 1.0E3
                if len(self.ws) == 1:
                    self.p_wt = self.p_wt[0]
            elif not self.version in self.versions:
                raise Exception("Version %s is not valid: version=[%s]"%(self.version, '|'.join(self.versions)))
        else:
            raise Exception("Version hasn't been set: version=[%s]"%('|'.join(self.versions)))
        return self.p_wt, self.u_wt, self.c_t
