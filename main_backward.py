import numpy as np

from src.ImpedancePy.shock_wave_compression.BackwardPropagationFromWindow import BackwardPropagationFromWindow
from src.ImpedancePy.shock_wave_compression.material_database.Kapton import Kapton
from src.ImpedancePy.shock_wave_compression.material_database.MgO import MgO
from src.ImpedancePy.shock_wave_compression.material_database.Quartz import Quartz

# Material.isentrope_calculator = ReflectedHugoniot()

measured_pressure = 315.03  # GPa
backward = BackwardPropagationFromWindow(materials=[Kapton(released=False),
                                                    MgO(is_stochastic=False, released=True),
                                                    Quartz(is_stochastic=False, released=True)])
backward.propagate(measured_pressure)
backward.plot()

ablator_intersections = np.array(
    [(intersection.particle_velocity, intersection.pressure) for intersection in
     backward.initial_intersections])
mean_ablator_pressure = np.mean(ablator_intersections[:, 1], dtype=np.float64)
std_ablator_pressure = np.std(ablator_intersections[:, 1], dtype=np.float64)
#
# print(f"Mean pressure: {mean_ablator_pressure}")
# print(f"Std pressure: {std_ablator_pressure}")
#
# plt.style.use('science')
# df = pd.DataFrame(ablator_intersections, columns=['Particle Velocity (km/s)', 'Pressure (GPa)'])
# sns.displot(data=df, x="Pressure (GPa)", kind="kde")
# plt.savefig(f'normal_2D_kde_{datetime.now()}.png')
# a=1
