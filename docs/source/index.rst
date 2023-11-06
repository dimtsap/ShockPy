.. Analytical shock-wave experiments documentation master file, created by
   sphinx-quickstart on Sun Sep  3 10:45:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ShockPy's documentation!
=============================================================

Table of contents
--------------------
* :ref:`General Info <General>`
* :ref:`Method <Method>`
* :ref:`Examples <Examples>`
* :ref:`Contents <Contents>`
* :ref:`Getting started <GettingStarted>`

.. _General:
General info
--------------------
This repository contains a library for modeling uncertainties in shock-wave experiments.
Specifically, using available experimental shock velocity and particle velocity measurements derived from shock-wave experiments,
Bayesian parameter estimation is performed to enhance the until now deterministic analytical Hugoniot formulas with uncertainty measures.

By combining multiple materials with uncertain Hugoniot equations and the Impedance Matching technique, shock-wave experiments can now be replicated using analytical formulas.
This approach provides experimentalists with statistical data about the possible experimental outcomes and thus expedites the experimentla deisgn process.


.. _Method:
Method
--------------------
Details of the methodology can be found in the `paper <TBD>`_.

.. _Examples:
Examples
--------------------
Inside the `documentation <https://dimtsap.github.io/ShockPy/auto_examples/index.html>`_ a set of examples can be found that illustrate the use of the code
for both forward and backward propagation of the shock-wave experiments allowing the option of including uncertainty to the predicted outcomes.
An illustration of forward experiment propagation for a three material experimental setup is provided below.

.. image:: _static/forward_experiment.png
   :scale: 20 %
   :alt: Forward propagation of shock-wave experiment using analytical equations
   :align: center


.. _Contents:
Contents
--------------------
The :code:`src` folder contains the implementation of the aforementioned methodology.
The are contained inside the :code:`material_database` module inside the :code:`src`.
There all the sample generated using the

.. _GettingStarted:
Getting started
--------------------
Users
^^^^^^^^^^^^^^
**1.** Create an Anaconda Python 3.9 virtual environment:

:code:`conda create -n shock_wave python==3.9`

:code:`conda activate shock_wave`


**2.** Install code and dependencies via the following commands:

:code:`pip install ShockPy`

Developers
^^^^^^^^^^^^^^
Given that the goal of this library is to integrate a wide variety of material with uncertain Hugoniot representations,
all contributions to the library are encouraged. To assist this process a :code:`.devcontainer` is include in the library that
allows for a replicable environment to be generated. This is environment is based on Centos 8 distribution due fortran
dependencies of one of the first-principles libraries integrated. Despite that the process is streamlined and requires no
further knowledge from the developers.

Developers are encouraged to use Visual Studio code as it offers a seamless integration with DevContainer as they are auto-detected and the user is directly prompted to open and work on them.
For further info on DevContainer you can refer `here <https://code.visualstudio.com/docs/devcontainers/containers>`_.


Mainteners
^^^^^^^^^^^^^^
Dimitris Tsapetis

:email: dtsapet1@jhu.edu


.. toctree::
   :hidden:

   Home <self>
   /material_database/index
   /material_states/index
   uncertainty_propagation
   bibliography.rst


