#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# try:
#     from setuptools import setup
# except ImportError:
#     from distutils.core import setup



#from setuptools import setup
#from setuptools import Extension
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
import os
import glob


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'sphinx',
    'sphinx-fortran',
    'numpy',
    'scipy',
    'pandas',
    'matplotlib',
    'PyYAML',
]

test_requirements = [
    'tox',
    'pytest',
    'coverall',
]

setup(
    name='fusedwake',
    version='0.1.0',
    description="A collection of wind farm flow models for FUSED-Wind",
    long_description=readme + '\n\n' + history,
    author="Pierre-Elouan Rethore",
    author_email='pire@dtu.dk',
    url='https://github.com/DTUWindEnergy/FUSED-Wake',
    packages=[
        'fusedwake',
    ],
    package_dir={'fusedwake':
                 'fusedwake'},
    include_package_data=True,
    install_requires=requirements,
    license="GNU Affero v3",
    zip_safe=False,
    keywords='fusedwake',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Affero v3',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    ext_modules=[Extension('fusedwake.gcl.fortran',
                           glob.glob(os.path.join('fusedwake', 'gcl', 'fortran',
                                                  'GCL.f'))),
                Extension('fusedwake.noj.fortran',
                           glob.glob(os.path.join('fusedwake', 'noj', 'fortran',
                                                  'NOJ.f'))),
                Extension('fusedwake.noj.mod_fortran',
                           glob.glob(os.path.join('fusedwake', 'noj', 'fortran',
                                                  'Mod_NOJ.f'))),
                Extension('fusedwake.bpa.fortran',
                           glob.glob(os.path.join('fusedwake', 'bpa', 'fortran',
                                                  'BPA.f')))],
)
