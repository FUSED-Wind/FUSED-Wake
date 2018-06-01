# import python.noj as noj
import fusedwake.noj.fortran as fnoj
import fusedwake.noj.fortran_mod as fnoj_mod
import numpy as np
from fusedwake.interface import BaseInterface


class NOJ(BaseInterface):
    # The different versions and their respective inputs
    inputs = {
        # 'py0': ['WF', 'WS', 'WD', 'K'],
        'fort_noj_av': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                        'av', 'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_noj': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                     'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_noj_s': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                       'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_mod_noj_av': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                            'av', 'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_mod_noj': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                         'rho', 'ws_ci', 'ws_co', 'ct_idle'],
        'fort_mod_noj_s': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'kj',
                           'rho', 'ws_ci', 'ws_co', 'ct_idle'],
    }
    # Default variables for running the wind farm flow model
    defaults = {
        'rho': 1.225,
        'K': 0.04,
        'version': 'fort_noj_s',
        'sup': 'quad',  # ['lin' | 'quad']
    }

    def fort_noj_av(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.kj = np.array([self.K])
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
            self.p_wt, self.t_wt, self.u_wt = fnoj.noj_av(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1:  # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fort_noj(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.kj = np.array([self.K])
        else:
            self.ws = self.WS
            self.wd = self.WD
            self.kj = np.full_like(self.wd, self.K)

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fnoj.noj(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1:  # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fort_noj_s(self):
        self.kj = self.K
        try:
            self.p_wt, self.t_wt, self.u_wt = fnoj.noj_s(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt

    def fort_mod_noj_av(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.kj = np.array([self.K])
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
            self.p_wt, self.t_wt, self.u_wt = fnoj_mod.mod_noj_av(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1:  # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fort_mod_noj(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.kj = np.array([self.K])

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fnoj_mod.mod_noj(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1:  # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fort_mod_noj_s(self):
        self.kj = self.K
        try:
            self.p_wt, self.t_wt, self.u_wt = fnoj_mod.mod_noj_s(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the followind inputs: {}, and the error message: {}'.format(
                self.version, self._get_kwargs(self.version), e))
        A = 0.25 * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt

    # def python_v0(self):
    #     self.p_wt, self.u_wt, self.c_t = noj.nojarsen_v0(**self._get_kwargs(self.version))
    #
    # def python_v1(self):
    #     self.p_wt, self.u_wt, self.c_t = noj.nojarsen(**self._get_kwargs(self.version))
