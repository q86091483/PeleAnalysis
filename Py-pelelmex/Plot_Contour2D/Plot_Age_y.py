#%%
import sys, os
import re
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib        as mpl

path_PeleAnalysis = os.path.abspath("../..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import imageio

# Input
# Where original plot files are stored
str_plane = "HRR_T_y=-1.000E-03"
data_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/"
              "Data_age/Slice2D_plt_lev0/conv3D/" + str_plane)
str_plane = "HRR_T_y=-1.000E-03"
der_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/"
              "Data_age/Slice2D_plt_lev0/conv3D/" + str_plane)
tfine = 5.960E-4
eps_age = 1E-6
# Case name
case_name = "conv3D"
plane_info = r"$\mathrm{Centre\; plane}$"
# Patterns of plotfiles to be processed
data_pattern = "plt_1*"
# Field_names
field_names = [["mixture_fraction", "mixture_fraction_userdef_0", "temp", "x_velocity"],
              ["age_0", "age_1", "mixAge_0", "mixAge_1"],
              ["agepv_0", "agepv_1", "mixAgepv_0", "mixAgepv_1"]]
#field_names = [["mixture_fraction", "temp"],["x_velocity", "y_velocity"]]
# Output data folder
fig_dir = os.path.abspath("/scratch/b/bsavard/zisen347/PeleAnalysis/Figure")
fig_slice_dir = os.path.join(fig_dir, "Slice2D_age_lev0")
# if print on screen or write to file
print_mode = 1
# Micromix case parameter
zst = 0.001
Djet = 4.5E-4
# Get plot file names
fnn = os.path.join(data_folder, data_pattern)
fns_unsorted = glob.glob(fnn)
def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".npz")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
nfns = len(fns_sorted)
# Get derived file names
fns_der_unsorted = glob.glob(os.path.join(der_folder, data_pattern))
def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".npz")[0]
  return float(sss)
fns_der_sorted = sorted(fns_der_unsorted, key=get_key)
# Plot parameter
npy = len(field_names); npx = len(field_names[0])
fig_unit_y = 1.5
labelsize = 14
loc_cb = [-1.0, 11.0, 5, 1.0]
loc_cb_long = [-1.0, 10.0, 7, 1.0]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# Colormap
ncmp = 128; ncmp_black = 16
arr_summer = mpl.cm.summer(np.arange(128))
arr_summer_black = np.ones((ncmp_black, 4))
arr_autumn = mpl.cm.autumn(np.arange(128))
arr_autumn_black = np.ones((ncmp_black, 4))
color_black = np.array([0.0416    , 0.        , 0.        , 1.])
for i in range(0, ncmp_black):
  r = float(i/ncmp_black)
  arr_summer_black[i,:] = r*arr_summer[ncmp_black,:] + (1-r)*color_black
  arr_autumn_black[i,:] = r*arr_autumn[ncmp_black,:] + (1-r)*color_black
arr_l = np.flipud(np.vstack((arr_autumn_black, arr_autumn)))
arr_Harbin = np.vstack((arr_l, arr_summer_black, arr_summer))
Harbin = mpl.colors.ListedColormap(arr_Harbin, name="Harbin", N = 2*(ncmp+ncmp_black))
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
for i, ifn in enumerate(range(0, nfns)):
  fn = fns_sorted[ifn]
  output_name = re.split("/", fn)[-1]
  output_name = re.split(".npz", output_name)[0]
  time = float(re.split("t=", output_name)[1])
  output_name = os.path.join(fig_plane_case_slice_dir, output_name)
  print(fns_sorted[ifn])
  print(fns_der_sorted[ifn])

  f = np.load(fns_sorted[ifn])
  fder = np.load(fns_der_sorted[ifn])
  x1D = f["x"]; Lx = x1D[-1] - x1D[0]
  y1D = f["y"]; Ly = y1D[-1] - y1D[0]
  extent=np.array([x1D[0], x1D[-1], 0, y1D[-1]]) / Djet
  mixf = np.zeros_like(f["density"])
  temperature = f["temp"]


  figsize = ((Lx/Ly)*fig_unit_y*npx*0.85, fig_unit_y*npy)
  fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)

  for ipy in range(0, npy):
    for ipx in range(0, npx):
      ax = axs[ipy, ipx]
      field_name = field_names[ipy][ipx]
      if (field_name == "HeatRelease"):
        vmin = 0.0; vmax = 1E11
        im=ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="hot",
                  extent=extent, aspect='equal')
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.text(x=-11, y=13, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black")
        if plane_info != None:
          ax.text(x=-11, y=15.5, s=plane_info, fontsize = labelsize-4)
        if time >= tfine:
          finest_level = 3
        else:
          finest_level = 2
        if True:
          ax.text(x=-11, y=18, s=r"$\mathrm{Finest \; level}="+str(finest_level)+"$", fontsize = labelsize-4)

        ax.set_title(r"$\mathrm{HRR}\;\mathrm{[J/m^3s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "HeatReleaseFI"):
        vmin = -1E11; vmax = 1E11
        im=ax.imshow(fder[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap=Harbin,
                  extent=extent, aspect='equal')
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent, linestyles=["-"])
        ax.text(x=12.0, y=11.0, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black")
        ax.set_title(r"$\mathrm{HRR}\;\mathrm{[J/m^3s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb_long, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, 0.0, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        #cb.ax.set_xticklabels(["%.1E"%vmax, "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar
        cb.ax.set_xticklabels([r"$+10^{11}$", r"$0$", r"+$10^{11}$"], color="white", fontsize=labelsize-6)  # horizontal colorbar
        ax.text(-2, 11.3, r"$\mathbf{FI<0}$", color="red", fontsize=6)
        ax.text(4, 11.3, r"$\mathbf{FI}>0}\;(\mathrm{Premixed})$", color="lime", fontsize=6)

      if (field_name == "temp"):
        vmin = 450; vmax = 2500
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[1E-3],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$T~\mathrm{[K]}$", fontsize=labelsize-2, pad=3)
        ax.text(x=12.0, y=11.0, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black")

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction"):
        vmin = 0; vmax = 1.1
        im = ax.imshow(fder[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z$", fontsize=labelsize-2, pad=3)
        ax.text(x=-11, y=13, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black", fontsize=16)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mag_vort"):
        vmin = 3.5; vmax = 6
        im = ax.imshow(np.log10(f[field_name]), origin="lower",
                  vmin = vmin, vmax = vmax, cmap="binary",
                  extent=extent, aspect='equal')
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.contour(mixf, levels=[0.12],
                   origin='lower', colors=['magenta'], extent=extent)

        ax.set_title(r"$\mathrm{log}(|\omega|)$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="red")
        cb.ax.yaxis.set_tick_params(color="red")
        cb.outline.set_edgecolor("red")
        cb.ax.set_xticklabels([str(vmin), str(vmax)], color="red", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "x_velocity"):
        vmin = -100; vmax = 100
        im = ax.imshow(f[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="seismic",
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_x \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_xticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "y_velocity"):
        vmin = -100; vmax = 100
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
        vmin = -250; vmax = 250
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
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
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
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

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


      if (ipy == npy-1):
        ax.set_xlabel(r"$x / D_\mathrm{j}$", fontsize = labelsize)
        ax.set_xticks(np.array([0, 5, 10, 15, 20]))
      else:
        #ax.set_xlabel([])
        ax.set_xticks(np.array([]))

      if (field_name == "mixture_fraction_userdef_0"):
        vmin = 0; vmax = 1.1
        im = ax.imshow(fder[field_name]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
#        ax.contour(mixf, levels=[zst],
#                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction_userdef_1"):
        vmin = 0; vmax = 1.1
        im = ax.imshow(fder[field_name]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
#        ax.contour(mixf, levels=[zst],
#                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z_\mathrm{air}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction_sum"):
        vmin = 0; vmax = 0.1
        phi = np.ones_like(fder["mixture_fraction_userdef_1"])
        phi = phi - fder["mixture_fraction_userdef_1"]/fder["density"]
        im = ax.imshow(phi-mixf, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
#        ax.contour(mixf, levels=[zst],
#                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z_\mathrm{air}+Z_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "age_0"):
        vmin = 0; vmax = 4.0E-4
        rhoMixf = fder["mixture_fraction_userdef_0"] + np.amax(fder["mixture_fraction_userdef_0"]) * eps_age
        mixf = fder["mixture_fraction_userdef_0"] / fder["density"]
        im = ax.imshow(fder[field_name] / rhoMixf, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[1-0.0623],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "age_1"):
        vmin = 0; vmax = 4.0E-4
        rhoMixf = fder["mixture_fraction_userdef_1"] + np.amax(fder["mixture_fraction_userdef_1"]) * eps_age
        mixf = fder["mixture_fraction_userdef_1"] / fder["density"]
        im = ax.imshow(fder[field_name] / rhoMixf, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')

        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{air}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "agepv_0"):
        vmin = 0; vmax = 1.0E-4
        rhoMixf = fder["mixture_fraction_userdef_0"] + np.amax(fder["mixture_fraction_userdef_0"]) * eps_age
        mixf = fder["mixture_fraction_userdef_0"] / fder["density"]
        im = ax.imshow(fder[field_name] / rhoMixf, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{p,jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "agepv_1"):
        vmin = 0; vmax = 1.0E-4
        rhoMixf = fder["mixture_fraction_userdef_1"] + np.amax(fder["mixture_fraction_userdef_1"]) * eps_age
        mixf = fder["mixture_fraction_userdef_1"] / fder["density"]
        im = ax.imshow(fder[field_name] / rhoMixf, origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        ax.contour(mixf, levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$\alpha_\mathrm{p,air}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixAge_0"):
        vmin = 0; vmax = 1.0E-3
        mixf = fder["mixture_fraction_userdef_0"] / f["density"]
        im = ax.imshow(fder["age_0"]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$Z_\mathrm{jet} \alpha_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixAge_1"):
        vmin = 0; vmax = 1.5E-3
        mixf = fder["mixture_fraction_userdef_1"] / f["density"]
        im = ax.imshow(fder["age_1"]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        print("Here: ", np.amax(mixf), zst)
        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$Z_\mathrm{air} \alpha_\mathrm{air}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixAgepv_0"):
        vmin = 0; vmax = 3.0E-4
        mixf = fder["mixture_fraction_userdef_0"] / f["density"]
        im = ax.imshow(fder["agepv_0"]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$Z_\mathrm{jet} \alpha_\mathrm{p,jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixAgepv_1"):
        vmin = 0; vmax = 3.0E-4
        mixf = fder["mixture_fraction_userdef_1"] / f["density"]
        im = ax.imshow(fder["agepv_1"]/fder["density"], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
        #ax.contour(mixf, levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.contour(temperature, levels=[1750],
                   origin='lower', colors=['magenta'], extent=extent)
        ax.set_title(r"$Z_\mathrm{air} \alpha_\mathrm{p,air}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar


      if (field_name == "forcing_age_0"):
        vmin = 0; vmax = 4E-6
        im = ax.imshow(fder[field_name], origin="lower",
                  vmin = vmin, vmax = vmax, cmap="jet",
                  extent=extent, aspect='equal')
#        ax.contour(mixf, levels=[zst],
#                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\mathrm{Forcing-} \alpha_\mathrm{jet}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_xticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar




      if (ipx == 0):
        ax.set_ylabel(r"$z / D_\mathrm{j}$", fontsize = labelsize)
        ax.set_yticks(np.array([0, 5, 10]))
      #else:
        #ax.set_xlabel([])
        #ax.set_yticks(np.array([]))

  if print_mode == 1:
    plt.savefig(output_name+".png", dpi=300, bbox_inches="tight")
    fig.clf()
    plt.close()
  print(output_name)
#%%
#with imageio.get_writer('/path/to/movie.gif', mode='I') as writer:
#    for filename in filenames:
#        image = imageio.imread(filename)
#        writer.append_data(image)