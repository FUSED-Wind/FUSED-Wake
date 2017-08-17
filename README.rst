===============================
FUSED-Wake
===============================

.. image:: https://img.shields.io/pypi/v/fusedwake.svg
        :target: https://pypi.python.org/pypi/fusedwake

.. image:: https://img.shields.io/travis/rethore/FUSED-Wake.svg
        :target: https://travis-ci.org/rethore/FUSED-Wake

.. image:: https://readthedocs.org/projects/fused-wake/badge/?version=master
        :target: https://fused-wake.readthedocs.org/en/latest/?badge=master
        :alt: Documentation Status


FUSED-Wake is a collection of wind farm flow models.

* Free software: GNU Affero v3 license
* Documentation: https://fused-wake.readthedocs.org/en/master/.

Notice
------
This package is in development, everything is in alpha mode. Expect dragons.

Features
--------
Currently FUSED-Wake has the following models implemented:

* G.C. Larsen model [Larsen_2009]_ (Python, Fortran)
* N.O. Jensen model [Jensen1983]_ (Fortran)
* Bastankhah & Porte-Agel model [Bastankhah_Porte-Agel_2014]_ (Fortran)

Roadmap
-------
The following models are planned to be added to this library:

* Stand-alone Dynamic Wake Meandering model [Keck_2015]_ (Python)
* Stationary G.C. Larsen [Larsen_2009]_ (Matlab)
* N.O. Jensen [NOJensen_1983]_ (Python, Matlab)
* Ainslie [Ainslie_1988]_ (Python, Fortran)
* Bastankhah & Porte-Agel model [Bastankhah_Porte-Agel_2014]_ (Python)
* Stand-alone Dynamic Wake Meandering model [Keck_2015]_ (Matlab)
* EllipSys3D Actuator Disk [Rethore_2013]_ (through REST-API)
* FUGA [Ott_2011]_ (through the Colonel module)

Dependencies
------------
This package has the following dependencies
* [windIO](https://github.com/rethore/windIO)
* numpy & scipy
* plotly (optional)
* jupyter (optional)
* pandas (optional)

Docs
----
Documentation is available online at https://fused-wake.readthedocs.org

You can build your own docs locally using the command

    $ make docs

Contribute
----------
See CONTRIBUTING_

Installation
------------

FUSED-Wake contains Fortran extensions that require a correctly configured Fortran compiler.

Windows compiler installation instructions
""""""""""""""""""""""""""""""""""""""""""

* Install Intel Fortran compiler, and activate as follows:

    $ "C:\\Program Files (x86)\\Intel\\Composer XE\\bin\\ifortvars.bat" intel64

or

* MinGW (instruction derived from `here <https://www.scivision.co/f2py-running-fortran-code-in-python-on-windows/>`_)

    1. Install numpy
    2. Install mingw-64 to ``c:\mingw with x86_64``, chose ``posix``, ``seh`` options
    3. Add MinGW bin folder (``C:\mingw\mingw64\bin``) to path variable
    4. Verify you can use gcc by typing gcc into Anaconda prompt
    5. Update your distutils configuration file to indicate you are using MinGW::

        [build]
        compiler=mingw32

    6. into either one of the following configuration files:

        * ``c:\Anaconda\Lib\distutils\distutils.cfg``
        * ``<user_folder>\AppData\Local\Continuum\Miniconda3\Lib\distutils\distutils.cfg``

Installing dependencies in a conda environment
""""""""""""""""""""""""""""""""""""""""""""""

Avoid conda and pip taking over packages from each other at random moments::

    conda install numpy scipy pandas jupyter plotly
    conda install -c conda-forge utm --no-deps
    pip install sphinx-fortran --no-deps

And the windIO dependency::

    git clone https://github.com/rethore/windIO.git
    cd windIO
    pip install -e ./ --no-deps

Finally, build and install FUSED-Wake::

    git clone https://github.com/DTUWindEnergy/FUSED-Wake.git
    cd FUSED-Wake
    pip install -e ./ --no-deps


Installing simply using pip
"""""""""""""""""""""""""""

::

    pip install numpy scipy pandas jupyter plotly utm sphinx-fortran

And the windIO dependency::

    git clone https://github.com/rethore/windIO.git
    cd windIO
    pip install -e ./

Finally, build and install FUSED-Wake::

    git clone https://github.com/DTUWindEnergy/FUSED-Wake.git
    cd FUSED-Wake
    pip install -e ./


Tests
-----
Local tests
"""""""""""
You can run the tests for your python environment using

  $ make tests

All tests
"""""""""
You can run all the tests for all the suported python versions

  $ make all-tests

Linting
"""""""
You can test if there are some flake8 issues

  $ make lint

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _CONTRIBUTING: _https://github.com/rethore/FUSED-Wake/blob/master/CONTRIBUTING.rst
