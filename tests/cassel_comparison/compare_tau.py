#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

angles = ["1p0", "3p5"]
cassel_names = ["cassel_tau_alpha{}.dat".format(a) for a in angles]
my_names = ["results_surf_alpha{}.dat".format(a) for a in angles]

cassel_data = list(map(lambda fn : np.genfromtxt(fn, delimiter=","), cassel_names))
my_data = list(map(lambda fn : np.genfromtxt(fn), my_names))

# Filtering first cassel data set
def filt(x, tau):
    for i in range(2, len(tau)-2):
        if x[i] <= -4.3 or x[i] > 10.0:
            tau[i] = (tau[i-2] + tau[i-1] + tau[i] + tau[i+1] + tau[i+2])/5.0

    return (x, tau)
cassel_data[0][:, 0], cassel_data[0][:, 1] = filt(cassel_data[0][:, 0], cassel_data[0][:, 1])

w = 4.0
h = w/1.33
fig1 = plt.figure(figsize=(w, h))
plt.plot(cassel_data[0][:,0], cassel_data[0][:,1], "ok", ms=5.0, mfc="None", mew=1.25)
plt.plot(my_data[0][:,0], my_data[0][:,1], "-k", lw=1.5)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()

fig2 = plt.figure(figsize=(w, h))
plt.plot(cassel_data[1][:,0], cassel_data[1][:,1], "ok", ms=5.0, mfc="None", mew=1.25, markevery=2, label="Cassel, et al.")
plt.plot(my_data[1][:,0], my_data[1][:,1], "-k", lw=1.5, label="Present Study")
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-15.0, 15.0])
plt.legend(loc=9, frameon=False)
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.tight_layout()

# Saving figures
path_dissertation = "/home/james/Documents/dissertation/document/images"
path_journal = "/home/james/Documents/mypapers/journal/incipient-sep/images"
fname1 = "comparison_cassel_1p0.pdf"
fname2 = "comparison_cassel_3p5.pdf"
fig1.savefig("{}/{}".format(path_dissertation, fname1))
fig1.savefig("{}/{}".format(path_journal, fname1))
fig2.savefig("{}/{}".format(path_dissertation, fname2))
fig2.savefig("{}/{}".format(path_journal, fname2))

#plt.show()
