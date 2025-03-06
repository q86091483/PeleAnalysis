#%%
import sys, os
import re
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt

path_PeleAnalysis = os.path.abspath("../..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import imageio

# Input
# Where original plot files are stored
#case_name = "GCI_jet3D"
case_name = "Ujet30_Re7200_rd53"
data_folder = ("/scratch/w47/zl5403/PeleAnalysis"
              "/Data/Slice2D_plt_lev0/" + case_name + "/HRR_T_y=0.000E+00")
# Case name
str_plane = "HRR_T_y=0.000E+00"
str_info = r"$\mathrm{Centreplane}$"
# Patterns of plotfiles to be processed
data_pattern = "plt_07125*"

# Field_names
field_names = [["Y(IC8H18)", "Y(NC7H16)", "Y(OH)"],
               ["Y(C7KET)", "Y(C8KET)", "Y(CH2O)"],
               ["x_velocity", "age_0", "mixture_fraction_userdef_0"],
               ["mixture_fraction", "age_1", "mixture_fraction_userdef_1"]]
field_names = [["Y(IC8H18)", "Y(NC7H16)", "mag_vort"],
               ["Y(C7KET)", "Y(C8KET)", "Y(CH2O)"],
               ["x_velocity", "y_velocity", "z_velocity"]]
field_names = [["Y(IC8H18)", "mixture_fraction", "mag_vort"],
               ["Y(C7KET)", "Y(C8KET)", "Y(CH2O)"],
               ["z_velocity", "y_velocity", "temp"],
               ["mixture_fraction_userdef_0", "age_0", "age_1"]]

#field_names = [["temp", "mag_vort"], ["Y(CH2O)", "Y(C7KET)"]]
# Output data folder
fig_dir = os.path.join(path_PeleAnalysis, "Figure")
fig_slice_dir = os.path.join(fig_dir, "Slice2D_plt_lev0")

#%%
print_mode = 1
zst = 0.0252
Djet = 1.7E-4
# Get file names
fns_unsorted = glob.glob(os.path.join(data_folder, data_pattern))
def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".npz")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
# Plot parameter
npy = len(field_names); npx = len(field_names[0])
fig_unit_y = 1.5
labelsize = 14
loc_cb = [4.0, 8.0, 12, 1.0]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#%%
if not os.path.exists(fig_dir):
  os.mkdir(fig_dir)
if not os.path.exists(fig_slice_dir):
  os.mkdir(fig_slice_dir)
fig_case_slice_dir = os.path.join(fig_slice_dir, case_name)
if not os.path.exists(fig_case_slice_dir):
  os.mkdir(fig_case_slice_dir)
fig_plane_case_slice_dir = os.path.join(fig_case_slice_dir, str_plane)
if not os.path.exists(fig_plane_case_slice_dir):
  os.mkdir(fig_plane_case_slice_dir)

#%%
for ifn, fn in enumerate(fns_sorted[0:]):
  output_name = re.split("/", fn)[-1]
  output_name = re.split(".npz", output_name)[0]
  time = float(re.split("t=", output_name)[1])
  output_name = os.path.join(fig_plane_case_slice_dir, output_name)

  f = np.load(fn)
  x1D = f["x"]; Lx = x1D[-1] - x1D[0]
  y1D = f["y"]; Ly = y1D[-1] - y1D[0]
  extent=np.array([x1D[0], x1D[-1], 0, y1D[-1]]) / Djet

  figsize = ((Lx/Ly)*fig_unit_y*npx*0.95, fig_unit_y*npy)
  fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)

  for ipy in range(0, npy):
    for ipx in range(0, npx):
      ax = axs[ipy, ipx]
      field_name = field_names[ipy][ipx]
      if (field_name == "HeatRelease"):
        vmin = 0.0; vmax = 5E11
        im=ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="hot",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.text(x=12.0, y=11.0, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="yellow")
        if str_info != None:
          ax.text(x=-7, y=14, s=str_info, fontsize = labelsize-6)
        ax.set_title(r"$\mathrm{HRR}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "temp"):
        vmin = 403; vmax = 2800
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$T~\mathrm{[K]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "mag_vort"):
        vmin = 0; vmax = 200000
        phi = np.log10(f[field_name])
        phi = f[field_name]
        im = ax.imshow(phi, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="binary",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\mathrm{log}(|\omega|)$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="red")
        cb.ax.yaxis.set_tick_params(color="red")
        cb.outline.set_edgecolor("red")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="red", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "x_velocity"):
        vmin = -30; vmax = 30
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="seismic",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_x \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), str(vmax)], color="white", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "y_velocity"):
        vmin = -30; vmax = 30
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="seismic",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_y \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "z_velocity"):
        vmin = -30; vmax = 30
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="seismic",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_z \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(NO)"):
        vmin = 0.0; vmax = 1E-4
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{NO}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(N2O)"):
        vmin = 0.0; vmax = 1E-5
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{N2O}$", fontsize=labelsize-2, pad=3)
        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(NNH)"):
        vmin = 0.0; vmax = 2E-7
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{NNH}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(N)"):
        vmin = 0.0; vmax = 2E-7
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{N}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(IC8H18)"):
        vmin = 0.0; vmax = 7E-1 * 0.3
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{i-C_8 H_{18}}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(NC7H16)"):
        vmin = 0.0; vmax = 3E-1 * 0.3
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{i-C_7 H_{16}}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(CH2O)"):
        vmin = 0.0; vmax = 1E-3
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="Greens",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{C H_{2} O}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(C7KET)"):
        vmin = 0.0; vmax = 5E-5
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="Blues",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{C7KET}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(C8KET)"):
        vmin = 0.0; vmax = 5E-5
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{C8KET}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar


      if (field_name == "Y(OH)"):
        vmin = 0.0; vmax = 2E-3
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_{\mathrm{OH}}}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction"):
        vmin = 0.0; vmax = 0.3
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction_userdef_0"):
        vmin = 0.0; vmax = 1.1
        im = ax.imshow(f[field_name] / f["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction_userdef_1"):
        vmin = 0.0; vmax = 1.1
        im = ax.imshow(f[field_name] / f["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z_\mathrm{cf}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "age_0"):
        vmin = 0.0; vmax = 5E-4
        rhoZ = f["mixture_fraction_userdef_0"]
        rhoZ = rhoZ + 1E-6 * np.amax(rhoZ)
        age = f[field_name] / rhoZ
        im = ax.imshow(f[field_name] / rhoZ, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "age_1"):
        vmin = 0.0; vmax = 3E-2
        rhoZ = f["mixture_fraction_userdef_1"]
        rhoZ = rhoZ + 1E-6 * np.amax(rhoZ)
        im = ax.imshow(f[field_name] / rhoZ, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{cf}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="black", fontsize=labelsize-6)  # horizontal colorbar



      if (ipy == npy-1):
        ax.set_xlabel(r"$x/D_{jet}$", fontsize = labelsize)
        #ax.set_xticks(np.array([0, 5, 10, 15, 20]))
      else:
        ax.set_xlabel(" ")
        ax.set_xticks(np.array([]))
  if print_mode == 1:
    plt.savefig(output_name+".png", dpi=500, bbox_inches="tight")
    fig.clf()
    plt.close()
  print(output_name)
#%%
#with imageio.get_writer('/path/to/movie.gif', mode='I') as writer:
#    for filename in filenames:
#        image = imageio.imread(filename)
#        writer.append_data(image)