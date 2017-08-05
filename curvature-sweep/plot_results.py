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
        self.r          = None
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
            line10 = lines[10].strip().split()
            self.alpha = float(line9[-1])
            self.r = float(line10[-1])

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
subdirs = ["r0p25", "r0p50", "r0p75", "r1p00"]

# Traversing directories and importing data
td_results = []
for d in dirs:
    for sd in subdirs:
        full_dir = "{}/{}".format(d, sd)
        td_results.append(td_result(full_dir))

# Filtering data
#td_plot_res = list(filter(lambda item : item.alpha%0.5 == 0.0, td_results))
td_plot_res = td_results

# Creating plot of scaled pressure
#clist = ["0.5" for i in range(0, 4)]
clist = ["#1f77b4" for i in range(0, 4)]
#clist = ["b", "g", "0.5"]
llist = ["-", "--", "-.", ":"]
w = 4.25
h = w/1.33
td_alpha_1 = list(filter(lambda item : item.alpha == 1.0, td_plot_res))
td_alpha_2 = list(filter(lambda item : item.alpha == 2.0, td_plot_res))
td_alpha_3 = list(filter(lambda item : item.alpha == 3.0, td_plot_res))
fig1a = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_1:
    plt.plot(p.x, p.p, llist[i], color=clist[0], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$p$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.legend()
plt.tight_layout()
fig1b = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_2:
    plt.plot(p.x, p.p, llist[i], color=clist[1], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$p$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()
fig1c = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_3:
    plt.plot(p.x, p.p, llist[i], color=clist[2], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
#plt.legend(loc=2, frameon=False)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$p$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()

# Creating plot of scaled shear stress
fig2a = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_1:
    plt.plot(p.x, p.tau, llist[i], color=clist[0], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()
fig2b = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_2:
    plt.plot(p.x, p.tau, llist[i], color=clist[1], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()
fig2c = plt.figure(figsize=(w, h))
plt.grid()
i = 0
for p in td_alpha_3:
    plt.plot(p.x, p.tau, llist[i], color=clist[2], label=r"$r\, = \, {:1.2f}$".format(p.r))
    i += 1
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\tau$", fontsize=14)
plt.xlim([-10.0, 10.0])
plt.tight_layout()

# Plotting the separation and reattachment points
td_sep_res = list(filter(lambda item : item.x_sep_rea != [], td_results))
td_alpha_2 = list(filter(lambda i : i.alpha == 2.0, td_sep_res))
td_alpha_3 = list(filter(lambda i : i.alpha == 3.0, td_sep_res))

fig3 = plt.figure(figsize=(w+1.0, h))
ax = fig3.add_subplot(111)
ax.plot([], "o", color="#1f77b4", mew=1.25, label=r"$x_S, \ \alpha=2$")
ax.plot([], "s", color="#1f77b4", mew=1.25, label=r"$x_R, \ \alpha=2$")
ax.plot([], "o", color="0.5", mew=1.25, label=r"$x_S, \ \alpha=3$")
ax.plot([], "s", color="0.5", mew=1.25, label=r"$x_R, \ \alpha=3$")
plt.grid()
for prof in td_alpha_2:
    ax.plot(prof.r, prof.x_sep_rea[0], "o", color="#1f77b4", mew=1.25)
    ax.plot(prof.r, prof.x_sep_rea[1], "s", color="#1f77b4", mew=1.25)
for prof in td_alpha_3:
    ax.plot(prof.r, prof.x_sep_rea[0], "o", color="0.5", mew=1.25)
    ax.plot(prof.r, prof.x_sep_rea[1], "s", color="0.5", mew=1.25)
plt.xlabel(r"$r$", fontsize=12)
plt.ylabel(r"$x$", fontsize=12)
plt.legend(loc='center right')
plt.tight_layout()

# Paths for images
path_dis = "/home/james/Documents/dissertation/document/images"
path_jou = "/home/james/Documents/mypapers/journal/incipient-sep/images"

# Saving images
fig1 = [fig1a, fig1b, fig1c]
fig2 = [fig2a, fig2b, fig2c]
fig1names = ["td_curv_a{:3d}_p.pdf".format(int(100.0*a)) for a in [1.0, 2.0, 3.0]]
fig2names = ["td_curv_a{:3d}_tau.pdf".format(int(100.0*a)) for a in [1.0, 2.0, 3.0]]
for f, fn in zip(fig1, fig1names):
    f.savefig("{}/{}".format(path_dis, fn))
    f.savefig("{}/{}".format(path_jou, fn))
for f, fn in zip(fig2, fig2names):
    f.savefig("{}/{}".format(path_dis, fn))
    f.savefig("{}/{}".format(path_jou, fn))
fig3name = "td_curv_sepsize.pdf"
fig3.savefig("{}/{}".format(path_dis, fig3name))
fig3.savefig("{}/{}".format(path_jou, fig3name))

# Closing figures
plt.close(fig1a)
plt.close(fig1b)
plt.close(fig1c)
plt.close(fig2a)
plt.close(fig2b)
plt.close(fig2c)

plt.show()
