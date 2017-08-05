#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import separationSize as sS

# Importing data
dirs = ["coarse", "medium", "fine"]
#mesh_dims = [r"$101\times\,51$", r"$201\times\,101$", r"$1001\times\,501$"]
mesh_dims = [r"$101\times\,51$", r"$201\times\,101$", r"$401\times\,201$"]
nnodes = [101*51, 201*101, 401*201]
d_list = [np.genfromtxt("{}/results_surf.dat".format(d)) for d in dirs]

# Finding the separation size for each data set
tau_zeros = [sS.getSepSize(d[:, 0], d[:, 1]) for d in d_list]
sep_size = [z[1] - z[0] for z in tau_zeros]

for d, s in zip(dirs, sep_size):
    pdiff = np.fabs(s - sep_size[0])/sep_size[0]*100.0
    print("{} \t {:1.6f} \t {:2.2f}".format(d, s, pdiff))

# Plotting wall shear stress
my_fs = 14
phi = 1.0 + 1.0/3.0
w = 5.25
h = w/phi
fig1 = plt.figure(figsize=(w, h))
lstyles = ["-", "--", "-."]
i = 0
for data in d_list:
    plt.plot(data[:, 0], data[:, 1], "{}k".format(lstyles[i]), alpha=0.575,
            label=mesh_dims[i], lw=1.5)
    i += 1
plt.xlabel(r"$x$", fontsize=my_fs)
plt.ylabel(r"$\tau_w$", fontsize=my_fs)
plt.tight_layout()

# Plotting pressure
fig2 = plt.figure(figsize=(w, h))
i = 0
for data in d_list:
    plt.plot(data[:, 0], data[:, 2], "{}k".format(lstyles[i]), alpha=0.575,
            label=mesh_dims[i], lw=1.5)
    i += 1
plt.xlabel(r"$x$", fontsize=my_fs)
plt.ylabel(r"$\tau_w$", fontsize=my_fs)
plt.legend(loc=2, fontsize=my_fs, frameon=False)
plt.tight_layout()

# Plotting separation size vs number of nodes
fig3 = plt.figure(figsize=(w, h))
rel_sep_size = [s/sep_size[0] for s in sep_size]
plt.plot(nnodes, rel_sep_size, "sw", mew=1.5)
plt.xlabel("Number of nodes", fontsize=my_fs)
plt.ylabel(r"$\frac{x_R - x_S}{(x_R - x_S)_{coarse}}$", fontsize=my_fs)
plt.xscale("log")
plt.tight_layout()

# Saving plots
jfm_dir = "/home/James/Desktop/myPapers/journal/jfm-incipient-separation/images"
#fig1.savefig("{}/td-mesh-study-tau.pdf".format(jfm_dir))
#fig2.savefig("{}/td-mesh-study-p.pdf".format(jfm_dir))

plt.show()
