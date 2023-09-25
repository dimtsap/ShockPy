Material Database
==================


This module contains functionality pertaining to the materials supported by this library:

All materials implementations are based on the baseclass :class:`.Material`.
The baseclass contains multiple methods that assist with the Impedance Matching procedure.
These methods include calculating the shocked stated of material for a given pressure, the isentrope initiating from a point on a material Hugoniot, as well as finding the intersection between Hugoniots and isentropes of consecutive materials.
A new material with either a deterministic or uncertain Hugoniot can be added in straightforward manner as analytically explained in the :doc:`User-defined material <new_materials>` page

The currently available materials are the following:

- :class:`.Kapton`: Class that produces Hugoniot of Kapton ablation material.

- :class:`.MgO`: Class that produces uncertain Hugoniot of MgO sample material.

- :class:`.Quartz`: Class that produces uncertain Hugoniot of Quartz window material.

.. toctree::
   :hidden:
   :maxdepth: 1

    Kapton <kapton>
    MgO <mgo>
    Quartz <quartz>
    User-defined material <new_materials>