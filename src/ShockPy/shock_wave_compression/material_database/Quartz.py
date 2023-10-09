import math

import numpy as np
from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
import os
import pickle


class Quartz(Material):

    def __init__(self, gamma_eff: float = 0.6,
                 initial_density: float = 2.65,
                 released=True,
                 is_stochastic=False):
        """
        Initializer of the Quartz class.

        :param gamma_eff: :math:`\Gamma_{eff}` is the updated Mie-Grüneisen (MG) parameter :math:`\Gamma`
         for the linear reference of the Mie-Grüneisen (MGLR) model, that assumes the release path can be accurately
         reproduced with constant :math:`\Gamma` along nearly its entirety. More details on the process of calculating
         the :math:`\Gamma_{eff}` given experimental data can be found in :cite:`Knudson`. Default value for MgO
         :math:`\Gamma_{eff}=0.6`.
        :param initial_density: Density of Quartz material at ambient environment conditions found
         `here <https://www.webmineral.com/data/Quartz.shtml>`_. Default Value :math:`2.65`.
         Units: :math:`g/cm^3`.
        :param is_stochastic: Boolean that defines if the Hugoniot used for this material is the deterministic one or
         multiple Hugoniots that are produced as a result of the Bayesian Inference framework will be used.
          * If :any:`True`, then samples of the Bayesian Inference will be used to generate the :code:`hugoniots_list` attributed.
          * If :any:`False`, only the nominal Hugoniot populates the list, and is derived from the least-square fitting of its analytical equation to the experimental data.
         Default value: :any:`False`.
        :param released: released: Boolean parameter indicating whether the reflected shock produced at the interface of the
         current material with the next is a rarefaction or a reshock.
          * If :any:`True`, then the algorithms assumes that the material undergoes release and compute the release part of the isentrope.
          * If :any:`False`, then the algorithm assumes the material will be reshocked and computed the respective part of the isentrope.
         Default value: :any:`True`.
        """
        super().__init__(gamma_eff, initial_density, is_stochastic, released)
        with open(os.path.join(os.path.dirname(__file__), "data", "samples_quartz_exponential0.p"), "rb") \
                as input_file:
            data = pickle.load(input_file)
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.278, 1.193, -2.505, -0.3701]))  # exponential hugoniot
        """
        This attribute represents the deterministic value all possible shocked states of the MgO material. 
        It is derived from the least-square fitting of its analytical equation to the experimental data.
        The fit of the data to a exponential curve yields the following coefficients. 
        :math:`\{a_0, a_1, a_2, a_3\} = \{6.278, 1.193, -2.505, -0.3701\}`
        """
        self.hugoniots_list = [self.nominal_hugoniot] if not is_stochastic \
            else [self.calculate_hugoniot(x) for x in data[:1000:10]]
        """
        A list of the possible material Hugoniots. If the initialization parameter :code:`is_stochastic` is False,
        then the nominal Hugoniot will be a member of the list. Otherwise, uncertain Hugoniots are 
        generated using samples produced by the Bayesian Inference and are contained in the 
        :code:`samples_quartz_exponential0.p` pickle file, inside the data folder of the materials_database module.
        """

    def analytical_shock_velocity_equation(self, parameters, hugoniot_particle_velocity):
        """
        Analytical Hugoniot equation for Quartz, given by the following exponential formula.
        :math:`U_s = a_0+a_1 \cdot u_p + a_2 \cdot u_p * e^{a_3 \cdot u_p}`

        :param parameters: A list containing the parameters :math:`[a_0, a_1, a_2, a_3]`.
        :param hugoniot_particle_velocity: The particle velocity :math:`(u_p)` domain where the values of the
         shock velocities :math:`U_s` of the MgO will be calculated.
        """
        return np.array([parameters[0] + parameters[1] * x + parameters[2] * x * math.exp(parameters[3] * x)
                         for x in hugoniot_particle_velocity])

