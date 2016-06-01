# import python.gau as gau
import fusedwake.gau.fortran as fgau
import numpy as np


class GAU(object):
    # The different versions and their respective inputs
    inputs = {
        # 'py0': ['WF', 'WS', 'WD', 'K'],
        'fort_gau_av': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ks','ng',
                  'av', 'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_gau': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ks','ng',
                  'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_gau_s': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ks','ng',
                  'rho', 'ws_ci', 'ws_co', 'ct_idle'],
    }
    # Default variables for running the wind farm flow model
    defaults = {
        'rho': 1.225,
        'K': 0.04,
        'version': 'fort_gau',
        'sup': 'quad', # ['lin' | 'quad']
        'NG': 4,
    }
    def __init__(self, **kwargs):
        self.set(self.defaults)
        self.set(kwargs)

    @property
    def versions(self):
        versions = list(self.inputs.keys())
        versions.sort()
        return versions

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
        if 'py' in version:
            return {k:getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}
        if 'fort' in version:
            # fortran only get lowercase inputs
            return {(k).lower():getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}

    def fort_gau_av(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.ks = np.array([self.K])
        if not hasattr(self, 'wt_available'):
            self.wt_available = np.ones([len(self.ws), self.WF.nWT])
            self.av = self.wt_available
        elif self.wt_available.shape == (len(self.ws), self.WF.nWT):
            self.av = self.wt_available
        else:
            # stacking up the availability vector for each flow case
            self.av = np.vstack([self.wt_available for i in range(len(self.ws))])

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fgau.gau_av(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                    self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1: # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fort_gau(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.ks = np.array([self.K])

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fgau.gau(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                    self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1: # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]


    def fort_gau_s(self):
        self.ks = self.K
        try:
            self.p_wt, self.t_wt, self.u_wt = fgau.gau_s(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                    self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt

    # def python_v0(self):
    #     self.p_wt, self.u_wt, self.c_t = gau.gauarsen_v0(**self._get_kwargs(self.version))
    #
    # def python_v1(self):
    #     self.p_wt, self.u_wt, self.c_t = gau.gauarsen(**self._get_kwargs(self.version))

    def __call__(self, **kwargs):
        self.set(kwargs)
        if hasattr(self, 'version'):
            getattr(self, self.version)()
            if not self.version in self.versions:
                raise Exception("Version %s is not valid: version=[%s]"%(self.version, '|'.join(self.versions)))
        else:
            raise Exception("Version hasn't been set: version=[%s]"%('|'.join(self.versions)))
        return self
