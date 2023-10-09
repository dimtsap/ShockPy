import os.path

import numpy as np
from ShockPy.shock_wave_compression.material_database.baseclass.Material import Material
import pickle


class MgO(Material):

    def __init__(self, gamma_eff: float = 1.1,
                 initial_density: float = 3.58,
                 is_stochastic: bool = False,
                 released: bool = True):
        """
        Initializer of the MgO class.

        :param gamma_eff: :math:`\Gamma_{eff}` is the updated Mie-Grüneisen (MG) parameter :math:`\Gamma`
         for the linear reference of the Mie-Grüneisen (MGLR) model, that assumes the release path can be accurately
         reproduced with constant :math:`\Gamma` along nearly its entirety. More details on the process of calculating
         the :math:`\Gamma_{eff}` given experimental data can be found in :cite:`Knudson`. Default value for MgO
         :math:`\Gamma_{eff}=1.1`.
        :param initial_density: Density of MgO material at ambient environment conditions found
         `here <https://www.americanelements.com/magnesium-oxide-1309-48-4>`_. Default Value :math:`3.58`.
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
        super().__init__(gamma_eff=gamma_eff,
                         initial_density=initial_density,
                         is_stochastic=is_stochastic,
                         released=released)
        with open(os.path.join(os.path.dirname(__file__), "data", "samplesMgO.p"), "rb") as input_file:
            data = pickle.load(input_file)
        self.nominal_hugoniot = self.calculate_hugoniot(np.array([6.6161, 1.4111, -0.016277]))
        """
        This attribute represents the deterministic value all possible shocked states of the MgO material. 
        It is derived from the least-square fitting of its analytical equation to the experimental data.
        The fit of the data to a quadratic polynomial yields the following coefficients. 
        :math:`\{a_0, a_1, a_2\} = \{6.6161, 1.4111, -0.016277\}`
        """
        self.hugoniots_list = [self.nominal_hugoniot] if not is_stochastic \
            else [self.calculate_hugoniot(x) for x in data[:1000:10]]
        """
        A list of the possible material Hugoniots. If the initialization parameter :code:`is_stochastic` is False,
        then the nominal Hugoniot will be a member of the list. Otherwise, uncertain Hugoniots are 
        generated using samples produced by the Bayesian Inference and are contained in the :code:`samplesMgO.p` 
        pickle file, inside the data folder of the materials_database module.
        """

    def analytical_shock_velocity_equation(self, parameters: list, hugoniot_particle_velocity: np.ndarray):
        """
        Analytical Hugoniot equation for MgO, given by the following second order polynomial.
        :math:`U_s = a_0+a_1 \cdot u_p + a_2 \cdot u_p^2`

        :param parameters: A list containing the parameters :math:`[a_0, a_1, a_2]`.
        :param hugoniot_particle_velocity: The particle velocity :math:`(u_p)` domain where the values of the
         shock velocities :math:`U_s` of the MgO will be calculated.
        """
        return parameters[0] + parameters[1] * hugoniot_particle_velocity + parameters[2] * np.square(
            hugoniot_particle_velocity)
