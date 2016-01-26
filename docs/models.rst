Models
======

Standalone Dynamic Wake Meandering
----------------------------------

The following is a general description of the standalone **Dynamic Wake Meandering (DWM)** implemented in FUSED Wake.
The present software is currently under development (beta version release).
It is an adaptation of the standalone Dynamic Wake Meandering model implementation of Rolf-Erik Keck ([Keck_2013b]_), based on the theory described in [Larsen_2008]_, [Larsen_2015]_, [Madsen_2009]_. 
As opposed to the implementation of the DWM model in the commercial aero-elastic software HAWC2 developed at DTU Wind Energy, the present version of the model is not coupled to an aero-elastic turbine model, and therefore, is only meant to be used in wind farm power production application.


Theory
^^^^^^
The Dynamic Wake Meandering is a widely used engineering model for simulating wind farm wake flows and the structural response of wind turbines in wind farm environment. Its computational demands are typically several order of magnitude lower than Computational Fluid Dynamics model. This model was first proposed by [Madsen_2003]_, with the ambition of combining an engineering flow model that captures the essential physics behind wake dynamics while maintaining low computational demand. This model therefore aims at performing wind farm design calculations, which typically involves a larger number of design iterations, while taking into account both load and power production aspects of each individual turbine. The combined load and power was, prior to the introduction of the DWM mode, not available in widely IEC standards models such as the wake model for load calculation of Frandsen [Frandsen_2007]_ and the industry standard Jensen's model [NOJensen_1983]_ for power production.
The original implementation of the model was proposed by [Madsen_2009]_ in DTU's aero-elastic code HAWC2 [Larsen_2007]_. The present version of the model is an extension of the standalone Matlab based model of
[Keck_2013b]_ applied to wind farm calculations where each turbine are modeled using a Blade Element Momemtum (BEM) approach.


.. image:: ../fusedwake/sdwm/DWM_Workflow.png


Text bkfaskfksdfds

.. math:: \sqrt{x^2}

* first item
* second item


Current development status
^^^^^^


fusedwake.sdwm package
======================

Submodules
----------

fusedwake.sdwm.DWM_GClarsenPicks module
---------------------------------------

.. automodule:: fusedwake.sdwm.DWM_GClarsenPicks
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.DWM_calc_mixL module
-----------------------------------

.. automodule:: fusedwake.sdwm.DWM_calc_mixL
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.DWM_flowfield_farm module
----------------------------------------

.. automodule:: fusedwake.sdwm.DWM_flowfield_farm
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.DWM_init_dict module
-----------------------------------

.. automodule:: fusedwake.sdwm.DWM_init_dict
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.DWM_main_BEM module
----------------------------------

.. automodule:: fusedwake.sdwm.DWM_main_BEM
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.DWM_misc module
------------------------------

.. automodule:: fusedwake.sdwm.DWM_misc
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.RUN_sDWM module
------------------------------

.. automodule:: fusedwake.sdwm.RUN_sDWM
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.cBEM module
--------------------------

.. automodule:: fusedwake.sdwm.cBEM
    :members:
    :undoc-members:
    :show-inheritance:

fusedwake.sdwm.cDWM module
--------------------------

.. automodule:: fusedwake.sdwm.cDWM
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: fusedwake.sdwm
    :members:
    :undoc-members:
    :show-inheritance:
