import math
import os

import numpy as np
import pandas as pd
current_dir = os.getcwd()


def velocities(params_array):
    number_of_samples = params_array.shape[0]
    results_list = []
    file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
    df = pd.read_excel(open(file_path, 'rb'), sheet_name="MgO_ordered")

    particle_velocity = df.Up.values.reshape(-1, 1)
    for i in range(number_of_samples):
        a0 = params_array[i][0]  # Initial Volume (A**3)
        a1 = params_array[i][1]  # Isothermal Bulk Modulus (GPa)
        a2 = params_array[i][2]  # Isothermal Pressure Derivative of Bulk Modulus (unitless)
        result = a0+a1*particle_velocity+a2*np.square(particle_velocity)
        results_list.append(result)

    return results_list


def velocities_quartz_exponential(params_array):
    number_of_samples = params_array.shape[0]
    results_list = []
    file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
    df = pd.read_excel(open(file_path, 'rb'), sheet_name="Sheet1")

    particle_velocity = df.Up.values.reshape(-1, 1)
    for i in range(number_of_samples):
        a0 = params_array[i][0]
        a1 = params_array[i][1]
        a2 = params_array[i][2]
        a3 = params_array[i][3]
        result = [a0 + a1 * x + a2 * x * math.exp(a3*x) for x in particle_velocity]
        results_list.append(result)

    return [[x[0] for x in results_list[0]]]


def velocities_quartz_cubic(params_array):
    number_of_samples = params_array.shape[0]
    results_list = []
    file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
    df = pd.read_excel(open(file_path, 'rb'), sheet_name="Sheet1")

    particle_velocity = df.Up.values.reshape(-1, 1)
    for i in range(number_of_samples):
        a0 = params_array[i][0]
        a1 = params_array[i][1]
        a2 = params_array[i][2]
        a3 = params_array[i][3]
        result = [a0 + a1 * x + a2 * x**2 + a3 * x ** 3 for x in particle_velocity]
        results_list.append(result)

    return [[x[0] for x in results_list[0]]]
