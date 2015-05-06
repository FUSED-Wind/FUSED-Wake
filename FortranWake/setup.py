#!/usr/bin/env python


from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import glob

setup(name='FortranWake',
      version='0.1',
      description='GCL and NOJ wind park wake models in Fortran,Python,FusedWind',
      author='Juan Pablo Murcia',
      author_email='jumu@dtu.dk',
      url='https://github.com/DTUWindEnergy/FUSED-Wake/FortranWake',
      keywords=['wind farm power', 'wind power plant power', 'wind turbine wakes',
                'wind farm flow model'],
      ext_modules=[Extension('FortranWake.GCL',
                             glob.glob(os.path.join('src',
                                                    'GCL.f'))),
                  Extension('FortranWake.NOJ',
                             glob.glob(os.path.join('src',
                                                    'NOJ.f')))],
      packages=['FortranWake'])
