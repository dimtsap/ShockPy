import os
from UQpy.Inference import *
from UQpy.RunModel import RunModel
from UQpy.SampleMethods import *
from UQpy.Inference import BayesParameterEstimation
from UQpy.Distributions import *
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

current_dir = os.getcwd()
file_path = os.path.join(current_dir, "MgO_hugoniot.xlsx")
df = pd.read_excel(open(file_path, 'rb'), sheet_name="MgO_ordered")

shock_velocity = df.Us.values.reshape(-1, 1)
particle_velocity = df.Up.values.reshape(-1, 1)
error_shock_Velocity = np.atleast_2d(df.errUs).reshape(-1, 1)
error_particle_Velocity = np.atleast_2d(df.errUp).reshape(-1, 1)

plt.figure()
plt.xlabel("Particle velocity (km/s)")
plt.ylabel("Shock velocity (km/s)")
plt.errorbar(particle_velocity.flatten(), shock_velocity.flatten(),
             xerr=error_particle_Velocity.flatten(),
             yerr=error_shock_Velocity.flatten(), fmt='o')

a0min, a0max = 4, 10
a1min, a1max = 0, 3
a2min, a2max = -0.25, 0.1

dist1 = Uniform(a0min, a0max - a0min)
dist2 = Uniform(a1min, a1max - a1min)
dist3 = Uniform(a2min, a2max - a2min)
#
dist_nominal = [dist1, dist2, dist3]
names_ = ['a0', 'a1', 'a2']

h_func = RunModel(model_script='pfn.py', model_object_name='velocities', var_names=names_)

prior = JointInd(marginals=[dist1, dist2, dist3])
inference_model = InferenceModel(nparams=3, runmodel_object=h_func, prior=prior)

print('Performing Bayesian parameter estimation..')

seed_x = np.zeros((10, 6))
x0 = np.array([6.6161, 1.4111, -0.016277])
sigma_0 = 0.15 * x0
for i5 in range(3):
    seed_x[:, i5] = x0[i5] * (1 + sigma_0[i5] * np.random.randn(10))

bayes_estimator = BayesParameterEstimation(data=shock_velocity, inference_model=inference_model, sampling_class=MMH,
                                           nsamples=10000, jump=20, nburn=20000, seed=x0)

s = bayes_estimator.sampler.samples
pickle.dump(s, open("src/shock_wave_compression/material_database/data/samplesMgO.p", "wb"))


import pickle
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity  # for the plots
import matplotlib
# matplotlib.use('TkAgg')
# plt.rc('text', usetex=True)
import statistics

s = pickle.load(open("src/shock_wave_compression/material_database/data/samplesMgO.p", "rb"))


# Function to plot posterior pdf from samples
def pdf_from_kde(domain, samples1d):
    bandwidth = 1.06 * np.std(samples1d) * samples1d.size ** (-1 / 5)
    kde = KernelDensity(bandwidth=bandwidth).fit(samples1d.reshape((-1, 1)))
    log_dens = kde.score_samples(domain)

    return np.exp(log_dens)


fig, ax = plt.subplots(1, 3, figsize=(20, 4))

domain = np.linspace(a0min, a0max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 0])
ax[0].plot(domain, pdf_, label='posterior')
ax[0].plot(domain, dist1.pdf(domain), label='prior')
ax[0].set_xlabel('$a_0$')
ax[0].set_ylabel('Probability Density')

domain = np.linspace(a1min, a1max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 1])
ax[1].plot(domain, pdf_, label='posterior')
ax[1].plot(domain, dist2.pdf(domain), label='prior')
ax[1].set_xlabel('$a_1$ (GPa)')

domain = np.linspace(a2min, a2max, 200)[:, np.newaxis]
pdf_ = pdf_from_kde(domain, s[:, 2])
ax[2].plot(domain, pdf_, label='posterior')
ax[2].plot(domain, dist3.pdf(domain), label='prior')
ax[2].set_xlabel("$a_2$")
plt.savefig("parameter_estimation.png")
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

plt.figure()
particle_velocities = [x[0] for x in np.linspace(0, 15, 200)[:, np.newaxis].tolist()]
shock_velocities = [m0 + m1 * x + m2 * x ** 2 for x in particle_velocities]

plt.plot(particle_velocities, shock_velocities)
plt.show()
