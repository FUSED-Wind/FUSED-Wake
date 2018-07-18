import numpy as np

class BaseInterface(object):
    def __init__(self, **kwargs):
        self.set(self.defaults)
        self.set(kwargs)

    @property
    def versions(self):
        versions = list(self.inputs.keys())
        versions.sort()
        return versions


    def update_position(self, pos):
        self.WF.update_position(pos)
        self.x_g, self.y_g, self.z_g = self.WF.get_T2T_gl_coord2()


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
            The version of the wind farm flow model to run
        """
        if 'py' in version:
            return {k:getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}
        if 'fort' in version:
            # fortran only get lowercase inputs
            return {(k).lower():getattr(self, k) for k in self.inputs[version] if hasattr(self, k)}

    def __call__(self, **kwargs):
        self.set(kwargs)
        if hasattr(self, 'version'):
            getattr(self, self.version)()
            if not self.version in self.versions:
                raise Exception("Version %s is not valid: version=[%s]"%(self.version, '|'.join(self.versions)))
        else:
            raise Exception("Version hasn't been set: version=[%s]"%('|'.join(self.versions)))
        return self
