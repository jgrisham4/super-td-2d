#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import separationSize as ss

class td_result:

    def __init__(self, directory):
        self.directory  = directory
        self.input_file = "input.inp"
        self.alpha      = None
        self.x          = None
        self.tau        = None
        self.p          = None
        self.x_sep_rea  = None  # Separation and reattachment points
        self.get_alpha()
        self.import_profiles()

    def get_alpha(self):
        if not os.path.isfile("{}/{}".format(self.directory, self.input_file)):
            print("\nWarning: Can't open {}/{}".format(self.directory, self.input_file))
        else:
            with open("{}/{}".format(self.directory, self.input_file), "r") as fp:
                lines = fp.readlines()
            line9 = lines[9].strip().split()
            self.alpha = float(line9[-1])

    def import_profiles(self):
        if not os.path.isfile("{}/results_surf.dat".format(self.directory)):
            print("\nWarning: Can't open {}/results_surf.dat.".format(self.directory))
        else:

            # Importing and separating data
            d        = np.genfromtxt("{}/results_surf.dat".format(self.directory))
            self.x   = d[:, 0]
            self.tau = d[:, 1]
            self.p   = d[:, 2]

            # Determining separation size
            self.x_sep_rea = ss.getSepSize(self.x, self.tau)


# Getting list of directories
dirs = sorted([f for f in os.listdir(os.getcwd()) if os.path.isdir(f)])

# Traversing directories and importing data
td_results = []
for d in dirs:
    td_results.append(td_result(d))

# Filtering data
td_plot_res = list(filter(lambda item : item.alpha%0.5 == 0.0, td_results))

# Creating plot of scaled pressure
clist = ["0.75" for i in range(0, 3)] + ["0.5" for i in range(0, 3)] + ["0.25"]
llist = ["-", "--", "-."] + ["-", "--", "-.", "-"]
w = 5.5
h = w/1.61
fig1 = plt.figure(figsize=(w, h))
i = 0
for p in td_plot_res:
    plt.plot(p.x, p.p, llist[i], color=clist[i], label=r"$\alpha = {:1.1f}$".format(p.alpha))
    i += 1
#plt.legend(loc=2, frameon=False)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$p$", fontsize=14)
plt.xlim([-20.0, 20.0])
plt.tight_layout()

# Creating plot of scaled shear stress
fig2 = plt.figure(figsize=(w, h))
i = 0
for p in td_plot_res:
    plt.plot(p.x, p.tau, llist[i], color=clist[i], label=r"$\alpha = {:1.1f}$".format(p.alpha))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-20.0, 20.0])
plt.legend(loc=3, frameon=False)
plt.tight_layout()

# Plotting the separation and reattachment points
td_sep_res = list(filter(lambda item : item.x_sep_rea != [], td_results))
fig3 = plt.figure(figsize=(w, h))
plt.plot([], "ok", mfc="None", mew=1.25, label=r"$x_S$")
plt.plot([], "sk", mfc="None", mew=1.25, label=r"$x_R$")
for prof in td_sep_res:
    plt.plot(prof.alpha, prof.x_sep_rea[0], "ok", mfc="None", mew=1.25)
    plt.plot(prof.alpha, prof.x_sep_rea[1], "sk", mfc="None", mew=1.25)
    plt.legend(loc=2, frameon=False)
plt.xlabel(r"$\alpha$", fontsize=14)
plt.ylabel(r"$x$", fontsize=14)
plt.tight_layout()

# Writing data to file
with open("td_results.dat", "w") as tdf:
    tdf.write("{:5} {:13} {:13}\n".format("alpha", "xS", "xR"))
    for obj in td_sep_res:
        tdf.write("{:1.3f} {:1.6e} {:1.6e}\n".format(obj.alpha, obj.x_sep_rea[0], obj.x_sep_rea[1]))

# Paths for images
path_dis = "/home/james/Documents/dissertation/document/images"
path_jou = "/home/james/Documents/mypapers/journal/incipient-sep/images"

# Saving images
fig1name = "cassel_p.pdf"
fig2name = "cassel_tau.pdf"
fig3name = "xS-triple-deck.pdf"
fig1.savefig("{}/{}".format(path_dis, fig1name))
fig1.savefig("{}/{}".format(path_jou, fig1name))
fig2.savefig("{}/{}".format(path_dis, fig2name))
fig2.savefig("{}/{}".format(path_jou, fig2name))
fig3.savefig("{}/{}".format(path_dis, fig3name))
fig3.savefig("{}/{}".format(path_jou, fig3name))

#plt.show()


