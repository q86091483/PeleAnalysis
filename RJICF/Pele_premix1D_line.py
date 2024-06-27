#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import cantera as ct
#gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/BurkeH2/chem.yaml")
gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/BurkeH2/chem.yaml")

# Input
# Where plt files are stored
plt_folder = "/scratch/b/bsavard/zisen347/scopingRuns/conv2D"
# Case name
case_name = "conv2D"
lref = 1.0
# Patterns of plotfiles to be processed
plt_pattern = "plt_01000*"
plane_x = np.array([])
plane_y = np.array([0.0]) * lref
plane_z = np.array([]) * lref
# Patterns of npz files to be plotted
plot_pattern = "plt_01000*"
plot_plane_x = np.array([])
plot_plane_y = np.array([0.0]) * lref
plot_plane_z = np.array([]) * lref
# Max level
max_level = 0
# Prefix
str_prefix = "HRR_T"
# Fields to be extracted
field_names = ["density", "temp", "mixture_fraction", "mag_vort",
               "x_velocity", "y_velocity", "z_velocity",
               "mixture_fraction_userdef_0", "age_0"]
for isp, spn in enumerate(gas_mix.species_names):
  field_names.append("Y("+spn+")")
# Fields to be plotted
plot_names = [["x_velocity", "Y(N2)"],
              ["age_0", "mixture_fraction_userdef_0"]]
# Output data folder
output_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Data_age"
fig_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Figure_age"
slice_name = "Slice2D_plt_lev"+str(max_level)
# Print mode
print_mode = 0
# Create Data folder
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
output_slice_dir = os.path.join(output_dir, slice_name)
if not os.path.exists(output_slice_dir):
  os.mkdir(output_slice_dir)
output_case_slice_dir = os.path.join(output_slice_dir, case_name)
if not os.path.exists(output_case_slice_dir):
  os.mkdir(output_case_slice_dir)
# Figure folder
if not os.path.exists(fig_dir):
  os.mkdir(fig_dir)
fig_slice_dir = os.path.join(fig_dir, slice_name)
if not os.path.exists(fig_slice_dir):
  os.mkdir(fig_slice_dir)
fig_case_slice_dir = os.path.join(fig_slice_dir, case_name)
if not os.path.exists(fig_case_slice_dir):
  os.mkdir(fig_case_slice_dir)

npy = len(plot_names); npx = len(plot_names[0])
fig_unit_y = 4.0
str_info = r"$\mathrm{Centreplane}$"
labelsize = 14
lw = 3.0
loc_cb = [0.0, 1E-3, 0.0, 1E-5]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#%%
# Input file names for writing npz data files
fns_unsorted = glob.glob(os.path.join(plt_folder, plt_pattern))
def get_key(s):
  ss = re.split("plt_", s)[-1]
  return int(ss)
fns_sorted = sorted(fns_unsorted, key=get_key)
for fn in fns_sorted:
  print("Going to read ", fn)

#%% Print npz data
for iy, pos in enumerate(plane_y):

  normal = 1

  folder_name = str_prefix + "_y=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):
    os.mkdir(folder_name)

  for ifn, fn in enumerate(fns_sorted):
    print("Processing iy =", iy, ", pos=", pos)

    normal = 1
    out_name = re.split("/", fn)[-1]
    time = HeaderData(fn).time
    out_name = out_name + "_t=" + "%.3E"%time
    out_name = os.path.join(folder_name, out_name)

    print(normal, fn, out_name)
    mand = Mandoline(fn,
                     fields=field_names,
                     limit_level=max_level,
                     serial=True,
                     verbose=1)
    mand.slice(normal=normal,
               pos=pos,
               outfile=out_name,
               fformat="array",
               uselog=1,
               )
#%%
# Plot
for iy, pos in enumerate(plane_y):

  normal = 1

  folder_name = str_prefix + "_y=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)

  fns_plot_unsorted = glob.glob(os.path.join(folder_name, plot_pattern))
  def get_plot_key(s):
    ss = re.split("t=", s)[-1]
    sss = re.split(".npz", ss)[0]
    return float(sss)
  fns_plot_sorted = sorted(fns_plot_unsorted, key=get_plot_key)

  str_plane = str_prefix+"="+"%.3E"%plane_y[iy]
  fig_plane_case_slice_dir = os.path.join(fig_case_slice_dir, str_plane)
  if not os.path.exists(fig_plane_case_slice_dir):
    os.mkdir(fig_plane_case_slice_dir)

  for ifn, fn in enumerate(fns_plot_sorted):
    #str_plane = str_prefix+"="+"%.3E"%plane_y[iy]
    #fig_plane_case_slice_dir = os.path.join(fig_case_slice_dir, str_plane)

    output_name = re.split("/", fn)[-1]
    output_name = re.split(".npz", output_name)[0]
    time = float(re.split("t=", output_name)[1])
    output_name = os.path.join(fig_plane_case_slice_dir, output_name)

    f = np.load(fn)
    x1D = f["x"]; Lx = x1D[-1] - x1D[0]
    y1D = f["y"]; Ly = y1D[-1] - y1D[0]
    extent=np.array([x1D[0], x1D[-1], y1D[0], y1D[-1]]) / lref

    figsize = (fig_unit_y*npx*1.25, fig_unit_y*npy)
    fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)

    for ipy in range(0, npy):
      for ipx in range(0, npx):
        ax = axs[ipy, ipx]
        field_name = plot_names[ipy][ipx]

        if (field_name == "density"):
          vmin = 1.0; vmax = 1.2
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$\rho \; \mathrm{[kg/m^3]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "x_velocity"):
          vmin = 0; vmax = 20
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$u_x \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "y_velocity"):
          vmin = -0.1; vmax = 0.1
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$u_y \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "z_velocity"):
          vmin = -0.1; vmax = 0.1
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$u_y \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "temp"):
          vmin = 299; vmax = 301
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$T \; \mathrm{[K]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(N2)"):
          vmin = 0.0; vmax = 1.0
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$Y_\mathrm{N2}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(O2)"):
          vmin = 0.0; vmax = 1.0
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$Y_\mathrm{O2}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(H2)"):
          vmin = 0.0; vmax = 1.0
          phi = f[field_name]
          phi_1D = np.mean(f[field_name], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$Y_\mathrm{H2}$", fontsize=labelsize-2, pad=3)

        if (field_name == "mixture_fraction_userdef_0"):
          vmin = 0.0; vmax = 1.2
          phi = f[field_name] / f["density"]
          phi_1D = np.mean(f[field_name] / f["density"], 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_ylim([vmin, vmax])
          ax.set_title(r"$mixf_1$", fontsize=labelsize-2, pad=3)

        if (field_name == "age_0"):
          vmin = 0.0; vmax = 8E-5
          phi = f[field_name] / f["mixture_fraction_userdef_0"]
          age_1D = phi / (f["age_0"] + 1E-6)
          phi_1D = np.mean(phi, 0)
          ax.plot(x1D, phi_1D, linewidth=lw, color="r")
          ax.set_xlim(np.array([vmin, vmax])*10)
          ax.set_ylim(np.array([vmin, vmax]))
          ax.set_title(r"$mixf_1$", fontsize=labelsize-2, pad=3)


        #if (ipy == npy-1):
          #ax.set_xlabel(r"$x$", fontsize = labelsize)
          #ax.set_xticks(np.array([0, 5, 10, 15, 20]))
        #else:
          #ax.set_xlabel([])
          #ax.set_xticks(np.array([]))
    if print_mode == 1:
      plt.savefig(output_name+".png", dpi=300, bbox_inches="tight")
      fig.clf()
      plt.close()
    print(output_name)

# %%
