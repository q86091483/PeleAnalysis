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
#gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/nuig_H2_4atm/chem.yaml")
gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/BurkeH2/chem.yaml")

# Input
# Where plt files are stored
plt_folder = "/scratch/b/bsavard/zisen347/scopingRuns/conv2D"
# Case name
case_name = "conv2D"
lref = 1.0
# Patterns of plotfiles to be processed
plt_pattern = "plt_00950*"
plane_x = np.array([])
plane_y = np.array([]) * lref
plane_z = np.array([0.0]) * lref
# Patterns of npz files to be plotted
plot_pattern = "plt_00950*"
plot_plane_x = np.array([])
plot_plane_y = np.array([]) * lref
plot_plane_z = np.array([0.0]) * lref
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
plot_names = [["x_velocity", "age_0"],
              ["mixture_fraction_userdef_0", "mixture_fraction"]]
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
fig_unit_y = 2.0
str_info = r"$\mathrm{Centreplane}$"
labelsize = 14
loc_cb = [0.0, 1E-3, 0.0, 1E-5]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plot_cbar = False
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
for iz, pos in enumerate(plane_z):

  normal = 2

  folder_name = str_prefix + "_z=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):
    os.mkdir(folder_name)

  for ifn, fn in enumerate(fns_sorted):
    print("Processing iz =", iz, ", pos=", pos)

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
for iz, pos in enumerate(plane_z):

  normal = 2

  folder_name = str_prefix + "_z=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)

  fns_plot_unsorted = glob.glob(os.path.join(folder_name, plot_pattern))
  def get_plot_key(s):
    ss = re.split("t=", s)[-1]
    sss = re.split(".npz", ss)[0]
    return float(sss)
  fns_plot_sorted = sorted(fns_plot_unsorted, key=get_plot_key)

  str_plane = str_prefix+"="+"%.3E"%plane_z[iz]
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

    figsize = ((Lx/Ly)*fig_unit_y*npx*0.95, fig_unit_y*npy*1.2)
    fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)

    for ipy in range(0, npy):
      for ipx in range(0, npx):
        ax = axs[ipy, ipx]
        field_name = plot_names[ipy][ipx]
        if (field_name == "HeatRelease"):
          vmin = 0.0; vmax = 1.2E11
          im=ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="hot",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.text(x=12.0, y=11.0, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="yellow")
          if str_info != None:
            ax.text(x=-7, y=14, s=str_info, fontsize = labelsize-6)
          ax.set_title(r"$\mathrm{HRR}$", fontsize=labelsize-2, pad=3)
          if plot_cbar:
            cax = ax.inset_axes(loc_cb, transform=ax.transData)
            cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                           ticks=[vmin, vmax])
            cb.ax.xaxis.set_tick_params(color="white")
            cb.ax.yaxis.set_tick_params(color="white")
            cb.outline.set_edgecolor("white")
            cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

        if (field_name == "mixture_fraction"):
          vmin = -0.2; vmax = 1.2
          im = ax.imshow(f[field_name], origin="lower",
                    vmin = vmin, vmax = vmax, cmap="jet",
                    extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$\mathrm{Postprocessed} \; Z_1$", fontsize=labelsize-2, pad=3)
          if plot_cbar:
            cax = ax.inset_axes(loc_cb, transform=ax.transData)
            cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                           ticks=[vmin, vmax])
            cb.ax.xaxis.set_tick_params(color="white")
            cb.ax.yaxis.set_tick_params(color="white")
            cb.outline.set_edgecolor("white")
            cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

        if (field_name == "temp"):
          vmin = 700; vmax = 2800
          im = ax.imshow(f[field_name], origin="lower",
                    vmin = vmin, vmax = vmax, cmap="jet",
                    extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$T~\mathrm{[K]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "mag_vort"):
          vmin = 3.5; vmax = 6
          im = ax.imshow(np.log10(f[field_name]), origin="lower",
                   vmin = vmin, vmax = vmax, cmap="binary",
                   extent=extent, aspect='equal')
          ax.contour(f["mixture_fraction"], levels=[zst],
                    origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$\mathrm{log}(|\omega|)$", fontsize=labelsize-2, pad=3)

        if (field_name == "x_velocity"):
          vmin = 8.0; vmax = 12
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="seismic",
                    extent=extent, aspect='equal')
          ax.set_title(r"$u_x \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

          #cax = ax.inset_axes(loc_cb, transform=ax.transData)
          #cb = fig.colorbar(im, cax=cax, orientation='horizontal',
          #                 ticks=[vmin, vmax])
          #cb.ax.xaxis.set_tick_params(color="white")
          #cb.ax.yaxis.set_tick_params(color="white")
          #cb.outline.set_edgecolor("white")
          #cb.ax.set_xticklabels([str(vmin), str(vmax)], color="white", fontsize=labelsize-6)  # horizontal colorbar
        if (field_name == "y_velocity"):
          vmin = -1; vmax =1
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="seismic",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          X, Y = np.meshgrid(x1D, y1D, )
          U = f["x_velocity"]
          V = f["y_velocity"]
          stream = ax.streamplot(X, Y, U, V,
              start_points = [(0.0001, 0.0001)],
                                  density=35)
          paths = stream.lines.get_paths()
          segments = stream.lines.get_segments()
          num_pts = len(segments)
          flow_line = np.full((num_pts, 2), np.nan)
          for i in range(num_pts):
            flow_line[i,:] = segments[i][0,:]
          ax.plot(flow_line[:,0], flow_line[:,1])
          #ax.scatter(flow_line[:,0], flow_line[:,1], c="green")
          ax.set_title(r"$u_y \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)
          #ax.set_xlim([0.0004, 0.00075])
          #ax.set_ylim([-0.0002, 0.0002])

        if (field_name == "z_velocity"):
          vmin = -1; vmax = 1
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="seismic",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$u_z \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(H2)"):
          vmin = 0.0; vmax = 0.07
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          #ax.set_title(r"$Y_\mathrm{H2}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(NO)"):
          vmin = 0.0; vmax = 1E-4
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$Y_\mathrm{NO}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(N2O)"):
          vmin = 0.0; vmax = 1E-5
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$Y_\mathrm{N2O}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(NNH)"):
          vmin = 0.0; vmax = 2E-7
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$Y_\mathrm{NNH}$", fontsize=labelsize-2, pad=3)

        if (field_name == "Y(N)"):
          vmin = 0.0; vmax = 2E-7
          im = ax.imshow(f[field_name], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[zst],
          #           origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$Y_\mathrm{N}$", fontsize=labelsize-2, pad=3)

        #if (ipy == npy-1):
          #ax.set_xlabel(r"$x$", fontsize = labelsize)
          #ax.set_xticks(np.array([0, 5, 10, 15, 20]))
        #else:
          #ax.set_xlabel([])
          #ax.set_xticks(np.array([]))
        if (field_name == "mixture_fraction_userdef_0"):
          vmin = -0.2; vmax = 1.2
          im = ax.imshow(f[field_name]/f["density"], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[0.0623],
          #          origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$\mathrm{Transported} \; Z_1$", fontsize=labelsize-2, pad=3)

        if (field_name == "age_0"):
          vmin = -0.0; vmax = 8.0E-5
          #im = ax.imshow(f[field_name]/f["mixture_fraction_userdef_0"], origin="lower",
          im = ax.imshow(f[field_name]/f["density"], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="jet",
                   extent=extent, aspect='equal')
          #ax.contour(f["mixture_fraction"], levels=[0.0623],
          #          origin='lower', colors=['white'], extent=extent)
          ax.set_title(r"$Z_1 \alpha_1 $", fontsize=labelsize-2, pad=3)


    if print_mode == 1:
      plt.savefig(output_name+".png", dpi=300, bbox_inches="tight")
      fig.clf()
      plt.close()
    print(output_name)

# %%
from scipy.interpolate import RegularGridInterpolator
fig, ax = plt.subplots()
X, Y = np.meshgrid(x1D, y1D, )
U = f["x_velocity"]
V = f["y_velocity"]
vmin = 0.0; vmax = 1.0
im = ax.imshow(f["mixture_fraction_userdef_0"]/f["density"], origin="lower",
                   vmin = vmin, vmax = vmax, cmap="binary",
                   extent=extent, aspect='equal')
stream = ax.streamplot(X, Y, U, V,
            start_points = [(0.0001, 0.0000)],
                           density=35)

paths = stream.lines.get_paths()
segments = stream.lines.get_segments()
num_pts = len(segments)
flow_line = np.full((num_pts, 2), np.nan)
for i in range(num_pts):
  flow_line[i,:] = segments[i][0,:]
  ax.plot(flow_line[:,0], flow_line[:,1], color="magenta", linewidth=3.0)

interp_age = RegularGridInterpolator((y1D, x1D), f["age_0"] / f["mixture_fraction_userdef_0"], "linear")
interp_u = RegularGridInterpolator((y1D, x1D), f["x_velocity"], "linear")
interp_v = RegularGridInterpolator((y1D, x1D), f["y_velocity"], "linear")

coords = np.stack([flow_line[:,1], flow_line[:,0]], axis=-1)
age_1D = interp_age(coords)
u_1D = interp_u(coords)
v_1D = interp_v(coords)
npt = u_1D.shape[0]

vel_1D = np.zeros_like(age_1D)
vel_1D = np.sqrt(u_1D**2 + v_1D**2)

dx_1D = np.zeros_like(age_1D)
dx_1D[0:-1] = np.sqrt((flow_line[1:,0] - flow_line[0:-1,0])**2 +
            (flow_line[1:,1] - flow_line[0:-1,1])**2)
dx_1D[-1] = 0.0

age_1D_uv = np.zeros_like(age_1D)
s_1D = np.zeros_like(age_1D)
for i in range(1, npt):
  age_1D_uv[i] = age_1D_uv[i-1] + dx_1D[i-1] / vel_1D[i-1]
  s_1D[i] = s_1D[i-1] + dx_1D[i-1]

lw = 2.0
fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(s_1D, age_1D, color="r", linewidth = lw, linestyle = "-",
        label = r"$\alpha_1$")
ax.plot(s_1D, age_1D_uv, color="b", linewidth = lw, linestyle = "--",
        label = r"$\int 1/|u| \mathrm{d}s$")
labelsize = 16
ax.set_xlabel(r"$\mathrm{Stream\; coordinate} \; s \; [m]$", fontsize = labelsize)
ax.set_ylabel(r"$\tau \; [s]$", fontsize = labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='both', which='minor', labelsize=labelsize)
ax.set_xlim([0, 9E-4])
ax.set_ylim([0, 9E-5])
ax.legend(fontsize = labelsize)


# %%
