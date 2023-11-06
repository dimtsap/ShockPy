Material States
==================

This module contains functionality pertaining to shocked and released material states observed during the course of a shock wave experiment.


The currently available material states are

- :class:`.Hugoniot`: Dataclass that contains a locus of points for all the possible shocked states of a material.

- :class:`.Isentrope`: Dataclass that contains the path along which a material is released from its shocked state. For the purposes of this library an isentropic release is considered to be a good enough approximation. More details on the calculation of the isentrope are provide by the two available implementations, namely :class:`.ReflectedHugoniot` and :class:`.IntegratedIsentrope`.

- :class:`.Intersection`: Dataclass that contains all the information of the intersection between an Isentrope and a Hugoniot.


.. toctree::
   :hidden:
   :maxdepth: 1

    Hugoniot <hugoniot>
    Isentrope <isentrope>