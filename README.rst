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
