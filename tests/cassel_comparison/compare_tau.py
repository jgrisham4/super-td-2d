#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt

angles = ["1p0", "3p5"]
cassel_names = ["cassel_tau_alpha{}.dat".format(a) for a in angles]
my_names = ["results_surf_alpha{}.dat".format(a) for a in angles]
my_names.append("results_surf_alpha3p5-old.dat")
my_names.append("results_surf_alpha3p5-new.dat")

cassel_data = list(map(lambda fn : np.genfromtxt(fn, delimiter=","), cassel_names))
my_data = list(map(lambda fn : np.genfromtxt(fn), my_names))

fig1 = plt.figure()
plt.plot(cassel_data[0][:,0], cassel_data[0][:,1], "ok", lw=1.5)
plt.plot(my_data[0][:,0], my_data[0][:,1], "--r", lw=1.5)
fig2 = plt.figure()
plt.plot(cassel_data[1][:,0], cassel_data[1][:,1], "ok", lw=1.5, label="Cassel's data")
plt.plot(my_data[1][:,0], my_data[1][:,1], "--g", lw=1.5, label="Initial")
plt.plot(my_data[2][:,0], my_data[2][:,1], "--b", lw=1.5, label="Restart 1")
plt.plot(my_data[3][:,0], my_data[3][:,1], "--r", lw=1.5, label="Restart 2")
plt.legend(loc=3)

plt.show()
