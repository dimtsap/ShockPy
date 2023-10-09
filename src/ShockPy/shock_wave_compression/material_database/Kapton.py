import os.path

import numpy as np

from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
from ShockPy.shock_wave_compression.material_states.Hugoniot import Hugoniot


class Kapton(Material):

    def __init__(self, gamma_eff: float = 0.6,
                 initial_density: float = 1.33,
                 released: bool = False):
        """
        Initializer of the Kapton class.

        :param gamma_eff: :math:`\Gamma_{eff}` is the updated Mie-Grüneisen (MG) parameter :math:`\Gamma`
         for the linear reference of the Mie-Grüneisen (MGLR) model, that assumes the release path can be accurately
         reproduced with constant :math:`\Gamma` along nearly its entirety. More details on the process of calculating
         the :math:`\Gamma_{eff}` given experimental data can be found in :cite:`Knudson`. Default value for Kapton
         :math:`\Gamma_{eff}=0.6`.
        :param initial_density: Density of Kapton material at ambient environment conditions found
         `here <https://dielectricmfg.com/knowledge-base/kapton/>`_. Default Value :math:`1.33`.
         Units: :math:`g/cm^3`.
        :param released: Boolean parameter indicating whether the reflected shock produced at the interface of the
         current material with the next is a rarefaction or a reshock.
          * If :any:`True`, then the algorithms assumes that the material undergoes release and compute the release part of the isentrope.
          * If :any:`False`, then the algorithm assumes the material will be reshocked and computed the respective part of the isentrope.
         Default value: :any:`False`.
        """
        super().__init__(gamma_eff, initial_density, released)
        hugoniot_P_up = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "kapton_hugoniot.txt"))
        
        p = hugoniot_P_up[:, 0]
        up = hugoniot_P_up[:, 1]
        us = p / (self.initial_density * up)
        v = self.initial_volume * (us - up) / us
        self.nominal_hugoniot = Hugoniot(pressures=p, particle_velocities=up, shock_velocities=us, volumes=v)
        """This attribute represents the deterministic value all possible shocked states of the Kapton material. 
        Since we did not have access to experimental data to learn the relationship between shock :math:`(U_s)` and
        particle velocities :math:`(u_p)` of the Kapton :class:`.Hugoniot`, it is extracted from the hugoniot command of 
        HYADES hydrocode :cite:`HYADES`. The extracted data can be found inside the :code:`kapton_hugoniot.txt`
        file inside the data folder of the material database module."""
        self.hugoniots_list: list[Hugoniot] = [self.nominal_hugoniot]
        """A list of possible material Hugoniots. Since Kapton is deterministic it contains the data for only one 
        :class:`.Hugoniot` curve."""

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        pass
