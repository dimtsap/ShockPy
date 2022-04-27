import os
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, DotProduct, RBF
import numpy as np

current_dir = os.getcwd()
#
# file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
# df = pd.read_excel(open(file_path, 'rb'), sheet_name="MgO_ordered")
#
# shock_velocity = df.Us.values.reshape(-1, 1)
# particle_velocity = df.Up.values.reshape(-1, 1)
# error_shock_Velocity = np.atleast_2d(df.errUs).reshape(-1, 1)
# error_particle_Velocity = np.atleast_2d(df.errUp).reshape(-1, 1)

# plt.figure()
# plt.xlabel("Particle velocity (km/s)")
# plt.ylabel("Shock velocity (km/s)")
# plt.errorbar(particle_velocity.flatten(), shock_velocity.flatten(),
#              xerr=error_particle_Velocity.flatten(),
#              yerr=error_shock_Velocity.flatten(), fmt='o')

from os import stat
import numpy as np
import matplotlib.pyplot as plt
from UQpy.Inference import *
from UQpy.RunModel import RunModel
from UQpy.SampleMethods import *
from scipy import interpolate
from UQpy.Inference import BayesParameterEstimation
from UQpy.Distributions import *
from sklearn.neighbors import KernelDensity  # for the plots
import pickle
import matplotlib

# a0min, a0max = 4, 10
# a1min, a1max = 0, 3
# a2min, a2max = -0.25, 0.1
#
# dist1 = Uniform(a0min, a0max - a0min)
# dist2 = Uniform(a1min, a1max - a1min)
# dist3 = Uniform(a2min, a2max - a2min)
# #
# dist_nominal = [dist1, dist2, dist3]
# names_ = ['a0', 'a1', 'a2']
#
# h_func = RunModel(model_script='pfn.py', model_object_name='velocities', var_names=names_)
#
# prior = JointInd(marginals=[dist1, dist2, dist3])
# inference_model = InferenceModel(nparams=3, runmodel_object=h_func, prior=prior)
#
# print('Performing Bayesian parameter estimation..')
#
# seed_x = np.zeros((10, 6))
# x0 = np.array([6.6161, 1.4111, -0.016277])
# sigma_0 = 0.15 * x0
# for i5 in range(3):
#     seed_x[:, i5] = x0[i5] * (1 + sigma_0[i5] * np.random.randn(10))
#
# bayes_estimator = BayesParameterEstimation(data=shock_velocity, inference_model=inference_model, sampling_class=MMH,
#                                            nsamples=10000, jump=20, nburn=20000, seed=x0)
#
# s = bayes_estimator.sampler.samples
# pickle.dump(s, open("dewaeleSamples0.p", "wb"))

file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
df = pd.read_excel(open(file_path, 'rb'), sheet_name="Sheet1")

particle_velocity_exp_Q = df.Up.values.reshape(-1, 1)
shock_velocity_exp_Q = df.Us.values.reshape(-1, 1)
pressure_exp_Q = df.Pressure.values.reshape(-1, 1)
density_exp_calc_Q = [3.584 * x[1] / (x[1] - x[0]) for x in zip(particle_velocity_exp_Q, shock_velocity_exp_Q)]


shock_velocity_exp_Q_list=[x[0] for x in shock_velocity_exp_Q]

b0min, b0max = -50, 100
b1min, b1max = -30, 100
b2min, b2max = -50, 50
b3min, b3max = -10, 10

dist4 = Uniform(b0min, b0max - b0min)
dist5 = Uniform(b1min, b1max - b1min)
dist6 = Uniform(b2min, b2max - b2min)
dist7 = Uniform(b3min, b3max - b3min)
#
dist_nominal = [dist4, dist5, dist6, dist7]
names_ = ['b0', 'b1', 'b2', 'b3']

# h_func = RunModel(model_script='pfn.py', model_object_name='velocities_quartz_exponential', var_names=names_)
h_func = RunModel(model_script='pfn.py', model_object_name='velocities_quartz_cubic', var_names=names_)

prior = JointInd(marginals=[dist4, dist5, dist6, dist7])
inference_model = InferenceModel(nparams=4, runmodel_object=h_func, prior=prior)

print('Performing Bayesian parameter estimation..')

seed_x = np.zeros((10, 4))
# x0 = np.array([6.278, 1.193, -2.505, -0.3701])
x0 = np.array([1.754, 1.862, -3-364e-02, 5.665e-4])
sigma_0 = 0.15 * x0
for i5 in range(4):
    seed_x[:, i5] = x0[i5] * (1 + sigma_0[i5] * np.random.randn(10))

bayes_estimator = BayesParameterEstimation(data=shock_velocity_exp_Q_list, inference_model=inference_model,
                                           sampling_class=MMH,
                                           nsamples=5000, jump=10, nburn=20000, seed=x0)

s = bayes_estimator.sampler.samples
pickle.dump(s, open("samples_quartz_cubic.p", "wb"))


from os import stat
import numpy as np
import matplotlib.pyplot as plt
from UQpy.Inference import *
from UQpy.RunModel import RunModel
from UQpy.SampleMethods import *
from scipy import interpolate
from UQpy.Inference import BayesParameterEstimation
from UQpy.Distributions import *
from sklearn.neighbors import KernelDensity  # for the plots
import pickle
import matplotlib

# import pickle
# import matplotlib.pyplot as plt
# import numpy as np
# from sklearn.neighbors import KernelDensity  # for the plots
# import matplotlib
# # matplotlib.use('TkAgg')
# # plt.rc('text', usetex=True)
# import statistics
#
# s = pickle.load(open("samplesMgO.p", "rb"))
#
#
# # Function to plot posterior pdf from samples
def pdf_from_kde(domain, samples1d):
    bandwidth = 1.06 * np.std(samples1d) * samples1d.size ** (-1 / 5)
    kde = KernelDensity(bandwidth=bandwidth).fit(samples1d.reshape((-1, 1)))
    log_dens = kde.score_samples(domain)

    return np.exp(log_dens)


# fig, ax = plt.subplots(1, 3, figsize=(20, 4))
#
# domain = np.linspace(a0min, a0max, 200)[:, np.newaxis]
# pdf_ = pdf_from_kde(domain, s[:, 0])
# ax[0].plot(domain, pdf_, label='posterior')
# ax[0].plot(domain, dist1.pdf(domain), label='prior')
# ax[0].set_xlabel('$a_0$')
# ax[0].set_ylabel('Probability Density')
#
# domain = np.linspace(a1min, a1max, 200)[:, np.newaxis]
# pdf_ = pdf_from_kde(domain, s[:, 1])
# ax[1].plot(domain, pdf_, label='posterior')
# ax[1].plot(domain, dist2.pdf(domain), label='prior')
# ax[1].set_xlabel('$a_1$ (GPa)')
#
# domain = np.linspace(a2min, a2max, 200)[:, np.newaxis]
# pdf_ = pdf_from_kde(domain, s[:, 2])
# ax[2].plot(domain, pdf_, label='posterior')
# ax[2].plot(domain, dist3.pdf(domain), label='prior')
# ax[2].set_xlabel("$a_2$")
# plt.savefig("parameter_estimation.png")
# plt.show()
#
# m0 = statistics.mean(s[:, 0])
# print(m0)
# s0 = statistics.stdev(s[:, 0])
# print(s0)
#
# m1 = statistics.mean(s[:, 1])
# print(m1)
# s1 = statistics.stdev(s[:, 1])
# print(s1)
#
# m2 = statistics.mean(s[:, 2])
# print(m2)
# s2 = statistics.stdev(s[:, 2])
# print(s2)
#
# plt.figure()
# particle_velocities = [x[0] for x in np.linspace(0, 15, 200)[:, np.newaxis].tolist()]
# shock_velocities = [m0 + m1 * x + m2 * x ** 2 for x in particle_velocities]
#
# plt.plot(particle_velocities, shock_velocities)
# plt.show()

#################################################
import pickle
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity  # for the plots
import matplotlib
# matplotlib.use('TkAgg')
# plt.rc('text', usetex=True)
import statistics

s = pickle.load(open("samples_quartz_cubic.p", "rb"))
a=1
plt.figure()
fig, ax = plt.subplots(1, 4, figsize=(20, 4))

domain = np.linspace(b0min, b0max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 0])
ax[0].plot(domain, pdf_, label='posterior')
ax[0].plot(domain, dist4.pdf(domain), label='prior')
ax[0].set_xlabel('$b_0$')
ax[0].set_ylabel('Probability Density')

domain = np.linspace(b1min, b1max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 1])
ax[1].plot(domain, pdf_, label='posterior')
ax[1].plot(domain, dist5.pdf(domain), label='prior')
ax[1].set_xlabel('$b_1$ (GPa)')

domain = np.linspace(b2min, b2max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 2])
ax[2].plot(domain, pdf_, label='posterior')
ax[2].plot(domain, dist6.pdf(domain), label='prior')
ax[2].set_xlabel("$b_2$")

domain = np.linspace(b3min, b3max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 3])
ax[3].plot(domain, pdf_, label='posterior')
ax[3].plot(domain, dist7.pdf(domain), label='prior')
ax[3].set_xlabel("$b_3$")

plt.savefig("parameter_estimation_quartz.png")
plt.show()

m0 = statistics.mean(s[:, 0])
print(m0)
s0 = statistics.stdev(s[:, 0])
print(s0)

m1 = statistics.mean(s[:, 1])
print(m1)
s1 = statistics.stdev(s[:, 1])
print(s1)

m2 = statistics.mean(s[:, 2])
print(m2)
s2 = statistics.stdev(s[:, 2])
print(s2)

m3 = statistics.mean(s[:, 3])
print(m3)
s3 = statistics.stdev(s[:, 3])
print(s3)

# plt.figure()
# particle_velocities = [x[0] for x in np.linspace(0, 15, 200)[:, np.newaxis].tolist()]
# shock_velocities = [m0 + m1 * x + m2 * x ** 2 for x in particle_velocities]
#
# plt.plot(particle_velocities, shock_velocities)
# plt.show()