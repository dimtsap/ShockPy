Isentrope
--------------------

The class Isentrope defined below is an dataclass that contains all relevant information for the release/reshock curve of a material.

Isentrope Class
^^^^^^^^^^^^^^^

Methods
~~~~~~~~~~~~~~~~~~
.. autoclass:: ImpedancePy.shock_wave_compression.material_states.Isentrope


Inside :py:mod:`ImpedancePy` there are two ways of calculating the release or reshock of one material, to reach the
impedance of the next during a shock-wave experiment. The first one is the rough approximation of the release isentrope
using the reflected Hugoniot. The second approximation is an integrated isentrope approach introduced for Quartz in
:cite:`Knudson`. The classes that enable both implementations are described below. Both classes can be considered as
static and provide two methods. One for calculation of the isentrope named :code:`calculate_isentrope` that aids in the
forward propagation of the shock-wave experiment, while the second  called :code:`find_previous_material_isentrope`
identifies the isentrope of the previous material that goes through the current material shocked state.



Reflected Hugoniot Class
^^^^^^^^^^^^^^^^^^^^^^^^

Methods
~~~~~~~~~~~~~~~~~~
.. autoclass:: ImpedancePy.shock_wave_compression.material_states.isentrope_calculators.ReflectedHugoniot
    :members: calculate_isentrope, find_previous_material_isentrope






Integrated Isentrope Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Methods
~~~~~~~~~~~~~~~~~~
.. autoclass:: ImpedancePy.shock_wave_compression.material_states.isentrope_calculators.IntegratedIsentrope
    :members: calculate_isentrope, find_previous_material_isentrope