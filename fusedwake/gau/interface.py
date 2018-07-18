# import python.gau as gau
import fusedwake.gau.fortran as fgau
import numpy as np
from fusedwake.interface import BaseInterface


class GAU(BaseInterface):
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
        else:
            self.ws = self.WS
            self.wd = self.WD
            self.ks = self.K

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
        self.ws = self.WS
        self.wd = self.WD
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


