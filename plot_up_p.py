import math
import pickle
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity  # for the plots
import matplotlib
# matplotlib.use('TkAgg')
# plt.rc('text', usetex=True)
import statistics
import os
import pandas as pd

current_dir = os.getcwd()

file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
df = pd.read_excel(open(file_path, 'rb'), sheet_name="MgO_ordered")

particle_velocity_exp = df.Up.values.reshape(-1, 1)
shock_velocity_exp = df.Us.values.reshape(-1, 1)
pressure_exp = df.Pressure.values.reshape(-1, 1)
density_exp_calc = [3.584 * x[1] / (x[1] - x[0]) for x in zip(particle_velocity_exp, shock_velocity_exp)]

s = pickle.load(open("samplesMgO.p", "rb"))

up = np.linspace(3, 20, 200)


for i in range(0, 1000, 5):
    a0 = s[i, 0]
    a1 = s[i, 1]
    a2 = s[i, 2]
    Us = [a0 + a1 * x + a2 * x * x for x in up]
    pressure = [3.584 * x[0] * x[1] for x in zip(up, Us)]
    plt.plot(up, pressure, color='gray', zorder=3)

Us_0 = [6.6161 + 1.4111 * x - 0.016277 * x * x for x in up]
pressure = [3.584 * x[0] * x[1] for x in zip(up, Us_0)]
plt.plot(up, pressure, color='green', zorder=2, label='MgO Hugoniot')
plt.scatter(particle_velocity_exp, pressure_exp, c="blue", zorder=3)

file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
df = pd.read_excel(open(file_path, 'rb'), sheet_name="Sheet1")

particle_velocity_exp_Q = df.Up.values.reshape(-1, 1)
shock_velocity_exp_Q = df.Us.values.reshape(-1, 1)
pressure_exp_Q = df.Pressure.values.reshape(-1, 1)
density_exp_calc_Q = [3.584 * x[1] / (x[1] - x[0]) for x in zip(particle_velocity_exp_Q, shock_velocity_exp_Q)]
s_Q = pickle.load(open("samples_quartz_exponential.p", "rb"))

for i in range(0, 1000, 5):
    b0 = s_Q[i, 0]
    b1 = s_Q[i, 1]
    b2 = s_Q[i, 2]
    b3 = s_Q[i, 3]
    Us = [b0 + b1 * x + b2 * x * math.exp(b3 * x) for x in up]
    pressure = [2.65 * x[0] * x[1] for x in zip(up, Us)]
    plt.plot(up, pressure, color='gray', zorder=1)

Us_0 = [1.754 + 1.862 * x - 3.364e-2 * x ** 2 + 5.666e-4 * x ** 3 for x in up]
pressure = [2.65 * x[0] * x[1] for x in zip(up, Us_0)]
plt.plot(up, pressure, color='magenta', zorder=2, label='Quartz Hugoniot')

# plt.scatter(particle_velocity_exp, pressure_exp, c="blue", zorder=3)
plt.scatter(particle_velocity_exp_Q, pressure_exp_Q, c="blue", zorder=3)

x = np.loadtxt('kapton_hugoniot.txt')
kapton_pressure = x[:131, 0]
kapton_up = x[:131, 1]
plt.plot(kapton_up, kapton_pressure, color='red', zorder=2, label='Kapton Hugoniot ')

plt.xlabel("Particle Velocity (g/cm^3)", fontsize=15)
plt.ylabel("Pressure (GPa)", fontsize=15)
plt.legend()
plt.show()