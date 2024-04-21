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
data_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/"
              "Py-pelelmex/Data/Slice2D/Micromix/HRR_T_z=2.025E-03")
# Case name
case_name = "Micromix"
str_plane = "HRR_T_z=2.025E-03" # Output folder nmae
str_info = r"$z/D_{jet}=4.5$"   # Reminder on the figure
# Patterns of plotfiles to be processed
data_pattern = "plt_*"
# Field_names
field_names = [["HeatRelease", "Y(NO)", "x_velocity"], 
               ["temp", "Y(N)", "y_velocity"],
               ["mag_vort", "Y(NNH)", "z_velocity"]]
#%% 
print_mode = 1 # print_mode=1 will close figure and savefig
zst = 0.0252
Djet = 4.5E-4
# Get file names
fns_unsorted = glob.glob(os.path.join(data_folder, data_pattern))
def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".npz")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
# Plot parameter
npy = len(field_names); npx = len(field_names[0])
fig_unit_y = 2.0
labelsize = 14
loc_cb = [-3.0, -2.0, 1.0, 4.0]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# Output data folder
fig_dir = os.path.abspath("../Figure")
#%%
if not os.path.exists(fig_dir):
  os.mkdir(fig_dir)
fig_slice_dir = os.path.join(fig_dir, "Slice2D")
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
        vmin = 0.0; vmax = 1.2E11
        im=ax.imshow(f[field_name], origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="hot", 
                  extent=extent, aspect='equal')
        #ax.contour(f["mixture_fraction"], levels=[zst],
        #           origin='lower', colors=['white'], extent=extent)
        ax.text(x=12.0, y=3.0, s=r"$t="+("%.3f"%(time*1000))+"\mathrm{[ms]}$", c="black")
        if str_info != None:
          ax.text(x=-3, y=5.0, s=str_info, fontsize = labelsize-6)
        ax.set_title(r"$\mathrm{HRR}$", fontsize=labelsize-2, pad=3)

        cax = ax.inset_axes(loc_cb, transform=ax.transData)
        cb = fig.colorbar(im, cax=cax, orientation='vertical',
                          ticks=[vmin, vmax])
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

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
        cb.ax.xaxis.set_tick_params(color="white")
        cb.ax.yaxis.set_tick_params(color="white")
        cb.outline.set_edgecolor("white")
        cb.ax.set_yticklabels([str(vmin), "%.1E"%vmax], color="white", fontsize=labelsize-6)  # horizontal colorbar

      if (field_name == "mag_vort"):
        vmin = 3.5; vmax = 6
        im = ax.imshow(np.log10(f[field_name]), origin="lower", 
                  vmin = vmin, vmax = vmax, cmap="binary", 
                  extent=extent, aspect='equal')
        ax.contour(f["mixture_fraction"], levels=[zst],
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
        vmin = -150; vmax = 150
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
        ax.set_xlabel(r"$x/D_{jet}$", fontsize = labelsize)
        ax.set_xticks(np.array([-0.0, 5, 10, 15, 20]))
      else:
        #ax.set_xlabel([])
        ax.set_xticks(np.array([]))
      if (ipx == 0):
        ax.set_ylabel(r"$y/D_{jet}$", fontsize = labelsize)
        ax.set_yticks(np.array([-4.0, 4.0]))
      else:
        #ax.set_xlabel([])
        ax.set_yticks(np.array([]))

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