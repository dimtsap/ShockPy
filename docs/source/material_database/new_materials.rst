.. _New_materials:
User-defined materials
-----------------------

All existing implementation of materials in this library are based on the abstract baseclass :class:`.Material`. This
class defines a number of supportive functions that are common between all material implementations.
These include mostly operations between the different material states, such as Hugoniots and isentropes and finding
their intersection to enable to Forward and backward propagation of the shock-wave experiment.

Material Class
^^^^^^^^^^^^^^

Methods
~~~~~~~~~~~~~~~~~~
.. autoclass:: ImpedancePy.shock_wave_compression.baseclass.Material
    :members: initialize_isentropes_at_pressure, calculate_intersection, find_previous_material_intersections_from_current

Attributes
~~~~~~~~~~~~~~~~~~
.. autoattribute:: ImpedancePy.shock_wave_compression.material_database.baseclass.Material.isentrope_calculator

Custom Material Example
^^^^^^^^^^^^^^^^^^^^^^^^

The code below provides an example of how to create a new material with a simple linear parametric Hugoniot equation.
The material in considered to be uncertain, as the distributions of its parameters are inferred using Bayesian Parameter Estimation.

>>> class CustomMaterial(Material):
>>>
>>>    def __init__(self,
>>>                 gamma_eff: float,
>>>                 initial_density:float,
>>>                 released:bool,
>>>                 is_stochastic:bool):
>>>        super().__init__(gamma_eff, initial_density, released)
>>>        with open(os.path.join(os.path.dirname(__file__), "data", "InferenceResult.p"), "rb") as input_file:
>>>            data = pickle.load(input_file)
>>>
>>>        self.nominal_hugoniot = self.calculate_hugoniot(np.array([1.3, 0.5]))
>>>
>>>        self.hugoniots_list = [self.nominal_hugoniot] if not is_stochastic \
>>>            else [self.calculate_hugoniot(x) for x in data[:1000:10]]
>>>
>>>
>>>    def analytical_shock_velocity_equation(self, parameters: list, hugoniot_particle_velocity: np.ndarray):
>>>        return parameters[0] + parameters[1] * hugoniot_particle_velocity