import python.gcl as gcl
import fortran as fgcl
import numpy as np


class GCL(object):
    """GC Larsen wake model applied to flat terrain/offshore wind farms [1].

    GCL(WF, WS, WD, TI, version, pars)

    This model is able to run with individual turbine rotor averaged WS,
    WD and TI. Different turbines types can be used in the same power
    plant; this includes different rotor diameter, operation, and heights.

    A modification of the rotor averaging quadrature was done. For more
    information see the jupyter notebook inside the test folder.

    This class wraps four different implementations (version):
    py_gcl_v0:  Python version. Double WT loop. Outter loop starts from
                the furthest downstream turbine up to the upstream ones.
                Inner loop estimates the wake due to each upstream turbine.
    py_gcl_v1:  Python version. Single WT loop from the further
                upstream turbines. The deficit each turbine produces in
                all the wake affected turbines is vectorized.
    fort_gcl:   Fortran version. Fastest implementation. It is able to
                receive multiple flowcases for fast computation of AEP.
                It can handle individual turbine availability (non
                operating turbines inside the plant) per flow case.
    fort_gclm:  Fortran version. It can handle individual undisturbed
                rotor averaged WS, WD, TI and multiple flowcases.
                It can handle individual turbine availability (non operating
                turbines inside the plant).

    Inputs
    ----------
    WF          fusedwake.WindFarm object
    WS          Individual rotor averaged undisturbed wind speed
    WD          Individual rotor averaged undisturbed wind direction
    TI          Individual rotor averaged ambient turbulent intensity
    version     py_gcl_v0, py_gcl_v1, fort_gcl, fort_gclm
    pars        GCL wake model parameters [a1, a2, a3, a4, b1, b2]

    Examples
    --------
    Definition of the workflow: wind farm and wake model
        gcl = GCL(WF=wf)

    For single flow case (Python versions can only run in this mode):
        gcl(WS=8, WD=270, TI=0.1)

    For multiple flowcases use array with same size for WS, WD and TI:
        WS_cases=np.arange(4,26)
        WD_cases=np.arange(0,360,2)
        WS_ms,WD_ms=np.meshgrid(WS_cases,WD_cases)
        WS=WS_ms.flatten()
        WD=WD_ms.flatten()
        gcl(WF=wf, WS=WS, WD=WD, TI=0.1*np.ones_like(WD), version='fort_gcl')

    For multiple flowcases and individual WS, WD and TI:
        WS_cases=np.arange(4,26)
        WD_cases=np.arange(0,360,2)
        WS_ms,WD_ms=np.meshgrid(WS_cases,WD_cases)
        WS=WS_ms.reshape(-1,1)*np.ones([1,wf.nWT])
        WS=WS+np.random.normal(0,0.5,size=WS.shape)
        WD=WD_ms.reshape(-1,1)*np.ones([1,wf.nWT])
        WD=WD+np.random.normal(0,3,size=WS.shape)
        gcl(WF=wf, WS=WS, WD=WD, TI=0.1*np.ones_like(WD), version='fort_gclm')

    @moduleauthor:: Juan P. Murcia <jumu@dtu.dk>
                    Pierre-Elouan Re'thore' <pire@dtu.dk>

    References:
    [1] Larsen GC. "A simple stationary semi-analytical wake model", 2009

    """
    # The different versions and their respective inputs
    inputs = {
        'py_gcl_v0': ['WF', 'WS', 'WD', 'TI', 'pars'],
        'py_gcl_v1': ['WF', 'WS', 'WD', 'TI', 'pars'],
        'fort_gcl': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ti',
                  'av', 'a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'rho', 'ws_ci',
                  'ws_co', 'ct_idle'],
        'fort_gclm': ['x_g', 'y_g', 'z_g', 'dt', 'p_c', 'ct_c', 'ws', 'wd', 'ti',
                  'av', 'a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'rho', 'ws_ci',
                  'ws_co', 'ct_idle'],
    }
    # Default variables for running the wind farm flow model
    defaults = {
        'rho': 1.225,
        'version': 'fort_gcl',
        'pars': [0.435449861, 0.797853685, -0.124807893, 0.136821858, 15.6298, 1.0],        'inflow': 'log',
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
        """ Set the attributes of a dictionary as instance variables.
        Prepares for the different versions of the wake model

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
            self.dt = np.array(self.WF.rotor_diameter)
            self.p_c = np.array(self.WF.power_curve)
            self.ct_c = np.array(self.WF.c_t_curve)
            self.ws_ci = np.array(self.WF.cut_in_wind_speed)
            self.ws_co = np.array(self.WF.cut_out_wind_speed)
            self.ct_idle = np.array(self.WF.c_t_idle)

    def _get_kwargs(self, version):
        """Prepare a dictionary of inputs to be passed to the wind farm flow model

        Parameters
        ----------
        version: str
            The version of the wind farm flow model to run ['py_gcl_v0' | 'py_gcl_v1' | 'fort0']
        """
        if 'py' in version:
            return {k:getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}
        if 'fort' in version:
            # fortran only get lowercase inputs
            return {(k).lower():getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}

    def fortran_gcl(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = np.array([self.WS])
            self.wd = np.array([self.WD])
            self.ti = np.array([self.TI])
        else:
            self.ws = self.WS
            self.wd = self.WD
            self.ti = self.TI
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

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fgcl.gcl_av(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the following inputs: {}, and the error message: {}'.format(
                    self.version, self._get_kwargs(self.version), e))
        A = 0.25 * np.pi * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1: # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def fortran_gclm(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.ws = self.WS*np.ones([1,self.WF.nWT])
            self.wd = self.WD*np.ones([1,self.WF.nWT])
            self.ti = self.TI*np.ones([1,self.WF.nWT])
        elif len(self.WD.shape) == 1:
            self.ws = self.WS.reshape([1,self.WF.nWT])
            self.wd = self.WD.reshape([1,self.WF.nWT])
            self.ti = self.TI.reshape([1,self.WF.nWT])
        else:
            self.ws = self.WS
            self.wd = self.WD
            self.ti = self.TI
        if not hasattr(self, 'wt_available'):
            self.wt_available = np.ones_like(self.ws)
            self.av = self.wt_available
        elif self.wt_available.shape == (self.ws.shape[0], self.WF.nWT):
            self.av = self.wt_available
        else:
            # stacking up the availability vector for each flow case
            self.av = np.vstack([self.wt_available for i in range(self.ws.shape[0])])
        self.a1 = self.pars[0] * np.ones(self.ws.shape[0])
        self.a2 = self.pars[1] * np.ones(self.ws.shape[0])
        self.a3 = self.pars[2] * np.ones(self.ws.shape[0])
        self.a4 = self.pars[3] * np.ones(self.ws.shape[0])
        self.b1 = self.pars[4] * np.ones(self.ws.shape[0])
        self.b2 = self.pars[5] * np.ones(self.ws.shape[0])

        # Run the fortran code
        try:
            self.p_wt, self.t_wt, self.u_wt = fgcl.gclm_av(**self._get_kwargs(self.version))
        except Exception as e:
            raise Exception('The fortran version {} failed with the following inputs: {}, and the error message: {}'.format(
                    self.version, self._get_kwargs(self.version), e))
        A = 0.25 * np.pi * self.WF.WT.rotor_diameter**2.0
        self.c_t = self.t_wt / (0.5 * A * self.rho * self.u_wt**2.0)
        self.p_wt *= 1.0E3  # Scaling the power back to Watt
        if len(self.ws) == 1: # We are only returning a 1D array
            self.p_wt = self.p_wt[0]
            self.u_wt = self.u_wt[0]
            self.c_t = self.c_t[0]

    def python_v0(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.WS = self.WS*np.ones([self.WF.nWT])
            self.WD = self.WD*np.ones([self.WF.nWT])
            self.TI = self.TI*np.ones([self.WF.nWT])
        else:
            self.WS = self.WS
            self.WD = self.WD
            self.TI = self.TI
        self.p_wt, self.u_wt, self.c_t = gcl.GCLarsen_v0(**self._get_kwargs(self.version))


    def python_v1(self):
        # Prepare the inputs
        if isinstance(self.WS, float) or isinstance(self.WS, int):
            self.WS = self.WS*np.ones([self.WF.nWT])
            self.WD = self.WD*np.ones([self.WF.nWT])
            self.TI = self.TI*np.ones([self.WF.nWT])
        else:
            self.WS = self.WS
            self.WD = self.WD
            self.TI = self.TI
        self.p_wt, self.u_wt, self.c_t = gcl.GCLarsen(**self._get_kwargs(self.version))

    def __call__(self, **kwargs):
        self.set(kwargs)
        if hasattr(self, 'version'):
            if   self.version == 'py_gcl_v0':
                self.python_v0()
            elif self.version == 'py_gcl_v1':
                self.python_v1()
            elif self.version == 'fort_gcl':
                self.fortran_gcl()
            elif self.version == 'fort_gclm':
                self.fortran_gclm()
            elif not self.version in self.versions:
                raise Exception("Version %s is not valid: version=[%s]"%(self.version, '|'.join(self.versions)))
        else:
            raise Exception("Version hasn't been set: version=[%s]"%('|'.join(self.versions)))
        return self
