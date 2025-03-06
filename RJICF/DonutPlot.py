#%%
import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = "./tmp"
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import cantera as ct
import pandas as pd

fig, axs = plt.subplots(1, 3, figsize=[8.5, 3])

ax = axs[0]
colors_gpu = ['limegreen', 'mediumturquoise', 'red', 'deepskyblue', 'yellow', "orange"]
labels_gpu = ["MAC Projection", "Scalar Advection", "Scalar Diffusion", "Scalar Reaction", "Velocity Advance", "Others"]
times_gpu = np.array([0.5, 0.103, 1.51, 4.36, 0.368, 0.83])
time_total_gpu = 0.0
wedges, txt1, txt2 = ax.pie(times_gpu, colors = colors_gpu, labels = [" "," "," "," "," "," "],
       autopct="%1.1f%%", pctdistance = 0.85,
       wedgeprops=dict(width=0.4))
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
  ang = (p.theta2 - p.theta1)/2. + p.theta1
  y = np.sin(np.deg2rad(ang))
  x = np.cos(np.deg2rad(ang))
  horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
  connectionstyle = f"angle,angleA=0,angleB={ang}"
  kw["arrowprops"].update({"connectionstyle": connectionstyle})
  #ax.annotate(labels_gpu[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
  #              horizontalalignment=horizontalalignment, **kw)
ax.set_title("4 CPU + 4 A100 GPU \n Narval")
ax.legend(labels_gpu, labels=labels_gpu, bbox_to_anchor=(-0.0, 1.05))
ax.text(-2.7, -1.4, r"Average time per step: " + "%.2f"%7.38 + " s", fontsize=14)

ax = axs[1]
colors_cpu = ['limegreen', 'mediumturquoise', 'red', 'deepskyblue', 'yellow', "orange"]
labels_cpu = ["MAC Projection", "Scalar Advection", "Scalar Diffusion", "Scalar Reaction", "Velocity Advance", "Others"]
times_cpu = np.array([1.2, 3.12, 15.04, 37.92, 2.27, 4.6])
time_total_cpu = 0.0
wedges, txt1, txt2 = ax.pie(times_cpu, colors = colors_cpu, labels = [" "," "," "," "," "," "],
       autopct="%1.1f%%", pctdistance = 0.85,
       wedgeprops=dict(width=0.4))
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
  ang = (p.theta2 - p.theta1)/2. + p.theta1
  y = np.sin(np.deg2rad(ang))
  x = np.cos(np.deg2rad(ang))
  horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
  connectionstyle = f"angle,angleA=0,angleB={ang}"
  kw["arrowprops"].update({"connectionstyle": connectionstyle})
  #ax.annotate(labels_gpu[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
  #              horizontalalignment=horizontalalignment, **kw)
ax.set_title("40 Intel Skylake CPUs \n Niagara")
ax.text(-0.2, -1.4,  r"%.1f"%77.53 + " s", fontsize=14)

ax = axs[2]
colors_cpu = ['limegreen', 'mediumturquoise', 'red', 'deepskyblue', 'yellow', "orange"]
labels_cpu = ["MACProjection()", "ScalarAdvection", "ScalarDiffusion", "ScalarReaction", "Velocity Advance", "Others"]
times_cpu = np.array([0.58, 2.3, 9.61, 22.08, 1.24, 4.8])
time_total_cpu = 0.0
wedges, txt1, txt2 = ax.pie(times_cpu, colors = colors_cpu, labels = [" "," "," "," "," "," "],
       autopct="%1.1f%%", pctdistance = 0.85,
       wedgeprops=dict(width=0.4))
bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
  ang = (p.theta2 - p.theta1)/2. + p.theta1
  y = np.sin(np.deg2rad(ang))
  x = np.cos(np.deg2rad(ang))
  horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
  connectionstyle = f"angle,angleA=0,angleB={ang}"
  kw["arrowprops"].update({"connectionstyle": connectionstyle})
  #ax.annotate(labels_gpu[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
  #              horizontalalignment=horizontalalignment, **kw)
ax.set_title("48 AMD EPYC CPUs \n Narval")
ax.text(-0.0, -1.4, r"%.1f"%45.01 + " s", fontsize=14)

plt.savefig("Donut.png", dpi=400, bbox_inches="tight")

#%%