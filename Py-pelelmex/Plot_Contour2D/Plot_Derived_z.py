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
str_plane = "HRR_T_z=2.025E-03"
data_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/"
              "Py-pelelmex/Data/Slice2D_plt_lev2/Micromix/" + str_plane)
str_plane = "HRR_T_z=2.025E-03"
der_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/"
              "Py-pelelmex/Data/Slice2D_der_lev2/Micromix/" + str_plane)
tfine = 5.960E-4
# Case name
case_name = "Micromix"
plane_info = r"$\mathrm{Horizotal\;plane}\;z/D_\mathrm{j}=4.5$"
str_plane = "HRR_T_z=2.025E-03" # Output folder nmae
# Patterns of plotfiles to be processed
data_pattern = "plt_*"
# Field_names
field_names = [["HeatRelease", "temp"], 
               ["HeatReleaseFI", "Y(NO)"],
               ["mixture_fraction", "y_velocity"]]
# Output data folder
fig_dir = os.path.abspath("../Figure")
fig_slice_dir = os.path.join(fig_dir, "Slice2D_der_lev2")
# if print on screen or write to file
print_mode = 1
# Micromix case parameter
zst = 0.0252
Djet = 4.5E-4
# Get plot file names
fns_unsorted = glob.glob(os.path.join(data_folder, data_pattern))
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
fig_unit_y = 2.0
labelsize = 22
loc_cb = [-3.5, -2.0, 1.0, 4.0]
loc_cb_long = [-3.0, -2.0, 3.0, 4.0]
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
for i, ifn in enumerate(range(0,nfns)):
  fn = fns_sorted[ifn]
  output_name = re.split("/", fn)[-1]
  output_name = re.split(".npz", output_name)[0]
  time = float(re.split("t=", output_name)[1])
  output_name = os.path.join(fig_plane_case_slice_dir, output_name)
  print(fns_sorted[ifn])
  print(fns_der_sorted[ifn])

  f = np.load(fn)
  fder = np.load(fns_der_sorted[ifn])
  y1D = f["y"]; Ly = y1D[-1] - y1D[0]
  x1D = f["x"]; Lx = x1D[-1] - x1D[0]
  extent=np.array([x1D[0], x1D[-1], y1D[0], y1D[-1]]) / Djet

  figsize = ((Lx/Ly)*fig_unit_y*npx*0.75, fig_unit_y*npy)
  fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)
  
  for ipy in range(0, npy):
    for ipx in range(0, npx):
      ax = axs[ipy, ipx]
      field_name = field_names[ipy][ipx]
      if (field_name == "HeatRelease"):
        vmin = 0.0; vmax = 1.0E11
        im=ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="hot", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.text(-10, y=5.5, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black", fontsize=labelsize-6)
        if plane_info != None:
          ax.text(x=-10, y=7.3, s=plane_info, fontsize = labelsize-6)
        if time >= tfine:
          finest_level = 3
        else: 
          finest_level = 2
        if True:
          ax.text(x=-10, y=9.1, s=r"$\mathrm{Finest \; level}="+str(finest_level)+"$", fontsize = labelsize-6)

        ax.set_title(r"$\mathrm{HRR}\;\mathrm{[J/m^3s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, 0, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([r"$+10^{11}$", r"$0$", r"+$10^{11}$"], color="white", fontsize=labelsize-6)  # horizontal colorbar


      if (field_name == "HeatReleaseFI"):
        vmin = -1E11; vmax = 1E11
        im = ax.imshow(fder[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap=Harbin, 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\mathrm{HRR}\;\mathrm{[J/m^3s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, 0, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([r"$+10^{11}$", r"$0$", r"+$10^{11}$"], color="white", fontsize=labelsize-8)  # horizontal colorbar
        ax.text(-3.5, -3.5, r"$\mathbf{FI<0}$", color="red", fontsize=labelsize-8)
        ax.text(-3.5, +3.0, r"$\mathbf{FI}>0}\;(\mathrm{Premixed})$", color="lime", fontsize=labelsize-8)

      if (field_name == "temp"):
        vmin = 700; vmax = 2800
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$T~\mathrm{[K]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        color_temp = "white"
        cb.ax.xaxis.set_tick_params(color=color_temp)
        cb.ax.yaxis.set_tick_params(color=color_temp)
        cb.outline.set_edgecolor(color_temp)
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color=color_temp, fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mixture_fraction"):
        vmin = 0.0; vmax = 1.0
        im = ax.imshow(fder[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        ax.contour(fder["mixture_fraction"], levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Z$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), str(vmax)], color="white", fontsize=labelsize-6)  # horizontal colorbar


      if (field_name == "mag_vort"):
        vmin = 3.5; vmax = 6
        im = ax.imshow(np.log10(f[field_name]), origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="binary", 
                  extent=extent, aspect='equal')
        ax.contour(fder["mixture_fraction"], levels=[zst],
                   origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$\mathrm{log}(|\omega|)$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="red")
        cb.ax.yaxis.set_tick_params(color="red")
        cb.outline.set_edgecolor("red")
        cb.ax.set_yticklabels([str(vmin), str(vmax)], color="red", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "x_velocity"):
        vmin = -250; vmax = 250
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="seismic", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst], 
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_x \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_yticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "y_velocity"):
        vmin = -50; vmax = 50
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="seismic", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst], 
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_y \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_yticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "z_velocity"):
        vmin = -250; vmax = 250
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="seismic", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst], 
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$u_z \; \mathrm{[m/s]}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="black")
        cb.ax.yaxis.set_tick_params(color="black")
        cb.outline.set_edgecolor("black")
        cb.ax.set_yticklabels([str(vmin), str(vmax)], color="black", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(NO)"):
        vmin = 0.0; vmax = 1E-4
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{NO}$", fontsize=labelsize-2, pad=3)

        loc_NO = [0,0,0,0]
        loc_NO[0] = loc_cb[0] - 15
        loc_NO[1] = loc_cb[1]
        loc_NO[2] = loc_cb[2]
        loc_NO[3] = loc_cb[3]
        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(N2O)"):
        vmin = 0.0; vmax = 1E-5
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{N2O}$", fontsize=labelsize-2, pad=3)
        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "Y(NNH)"):
        vmin = 0.0; vmax = 2E-7
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst], 
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{NNH}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar
      if (field_name == "Y(N)"):
        vmin = 0.0; vmax = 2E-7
        im = ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="jet", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst], 
        #           origin='lower', colors=['white'], extent=extent)
        ax.set_title(r"$Y_\mathrm{N}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar


      if (ipy == npy-1):
        ax.set_xlabel(r"$x/D_{j}$", fontsize = labelsize)
        ax.set_xticks(np.array([-0.0, 5, 10, 15, 20]))
      else:
        #ax.set_xlabel([])
        ax.set_xticks(np.array([]))
      if (ipx == 0):
        ax.set_ylabel(r"$y/D_{j}$", fontsize = labelsize)
        ax.set_yticks(np.array([-4.0, 4.0]))
      else:
        #ax.set_xlabel([])
        ax.set_yticks(np.array([]))

      ax.tick_params(axis='both', which='major', labelsize=labelsize)
      ax.tick_params(axis='both', which='minor', labelsize=labelsize)
  #fig.tight_layout()
  plt.subplots_adjust(wspace=0.05, hspace=0.01)
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