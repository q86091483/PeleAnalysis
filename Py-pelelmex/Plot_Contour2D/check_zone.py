#%%
import sys, os
import re
import matplotlib
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

slice_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/Slice2D_derived/Micromix"
fn_y = os.path.join(slice_folder, "Derived_y=-9.000E-04/plt_11800_derived_t=6.213E-04.npz")
fn_x0 = os.path.join(slice_folder, "Derived_x=9.000E-04/plt_11800_derived_t=6.213E-04.npz")
fn_x1 = os.path.join(slice_folder, "Derived_x=2.250E-03/plt_11800_derived_t=6.213E-04.npz")
fn_x2 = os.path.join(slice_folder, "Derived_x=3.600E-03/plt_11800_derived_t=6.213E-04.npz")
fn_x3 = os.path.join(slice_folder, "Derived_x=6.750E-03/plt_11800_derived_t=6.213E-04.npz")
fn_x4 = os.path.join(slice_folder, "Derived_x=1.080E-02/plt_11800_derived_t=6.213E-04.npz")

xps = np.array([9E-4, 2.25E-3, 6.75E-3, 1.08E-2])

xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3
Djet = 4.5E-04
#%%
extent=np.array([xmin, xmax, zmin, zmax]) / Djet
extent_x = np.array([ymin, ymax, zmin, zmax]) / Djet 
tet = np.load(fn_y)
hrr = tet["HeatRelease"]
mixfrac = tet["mixture_fraction"]
pv = tet["pv"]
rho = tet["rho"]
FI = tet["FI"]

tet_x0 = np.load(fn_x0)
tet_x1 = np.load(fn_x1)
tet_x2 = np.load(fn_x2)
tet_x3 = np.load(fn_x3)
tet_x4 = np.load(fn_x4)
tet_xs = [tet_x0, tet_x1, tet_x2, tet_x3]

print("Max/min rho: ",  np.amin(rho), np.amax(rho))
print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min progvar: ",  np.amin(pv), np.amax(pv))
print("Max/min FI: ",  np.amin(FI), np.amax(FI))

#
A = 0.55
B = 0.2
J = 6.6
Djet = 4.5E-4
xs = np.linspace(0, 0.01, 1000)
ys = J * Djet * A * np.power(xs / (J*Djet), B)
def get_trajectory_z(x):
  return J * Djet * A * np.power(x / (J*Djet), B) 
#%%
# Z - y
fig, ax = plt.subplots()
vmin = 0; vmax = 1.0
im = ax.imshow(mixfrac, origin="lower", 
                vmin=vmin, vmax=vmax, 
                cmap="viridis", extent=extent)
ax.plot(xs/Djet, ys/Djet, color="magenta", linewidth=3.0)
ax.set_xlabel(r"$x/D_j$", fontsize=16)
ax.set_ylabel(r"$z/D_j$", fontsize=16)
ax.set_title(r"$\mathrm{Mixture \; fraction}$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig("z_trajectory.png", dpi=300, bbox_inches="tight")
#%%
# hrr - y
fig, ax = plt.subplots()
vmin = 0; vmax = 1.0E11
im = ax.imshow(hrr, origin="lower", 
                vmin=vmin, vmax=vmax, 
                cmap="hot", extent=extent)
ax.plot(xs/Djet, ys/Djet)
divider = make_axes_locatable(ax)
ax.contour(tet["mixture_fraction"], levels=[0.0252],
          origin='lower', 
          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
zone = tet_xs[3]["zone"]
vmin = -1.8E-3; vmax = 1.8E-3
ax.imshow(zone, origin="lower")
im = ax.imshow(zone, origin="lower", aspect='equal',
                vmin=vmin, vmax=vmax, 
                cmap="jet", extent=extent_x)
#%%
ncols = 2; nrows = 2
field_name = "HeatRelease"
vmin = 0.0; vmax = 1E11
cmap = "hot"
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(4, 6))
k0 = 1.5
xl = np.linspace(-4, 0, 50)
xr = np.linspace( 0, 4, 50)
plot_slope=False
for ipx in range(0, ncols):
  for ipy in range(0, nrows):
    ix = ipx * ncols + ipy
    ax = axs[ipx,ipy]
    phi = tet_xs[ix][field_name]
    mf = tet_xs[ix]["mixture_fraction"]
    im = ax.imshow(phi, origin="lower", aspect='equal',
                vmin=vmin, vmax=vmax, 
                cmap=cmap, extent=extent_x)
    ax.contour(mf, levels=[0.0252], aspect='equal',
          origin='lower', 
          colors=['cyan'], extent=extent_x)
    xloc = xps[ix]
    zloc = get_trajectory_z(xloc)

    y0 = -2 
    z0 = zloc / Djet
    b0 = z0 - k0 * y0 
    if plot_slope:
      ax.plot(xl, k0*xl+b0,color="cyan", linestyle="--")
    y0 = -2 
    z0 = zloc / Djet
    b0 = z0 - (-k0) * y0 
    if plot_slope:
      ax.plot(xl, (-k0)*xl+b0,color="cyan", linestyle="--")
    y0 = 2 
    z0 = zloc / Djet
    b0 = z0 - k0 * y0 
    if plot_slope:
      ax.plot(xr, k0*xr+b0,color="cyan", linestyle="--")
    y0 = 2 
    z0 = zloc / Djet
    b0 = z0 - (-k0) * y0 
    if plot_slope:
      ax.plot(xr, (-k0)*xr+b0,color="cyan", linestyle="--")
    ax.set_ylim([0, 12])
    ax.scatter(-2, zloc/Djet, color="magenta")
    ax.scatter(+2, zloc/Djet, color="magenta")
    ax.set_title(r"$x/D_j=" + str(int(xloc/Djet))+"$")



    if False:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes('right', size='5%', pad=0.05)
      fig.colorbar(im, cax=cax, orientation='vertical')
    if (ipx == nrows-1):
      ax.set_xlabel(r"$z / D_\mathrm{j}$", fontsize = 14)
      ax.set_xticks(np.array([-4, -2, 0, 2, 4]))
    else:
        #ax.set_xlabel([])
      ax.set_xticks(np.array([]))

    if (ipy == 0):
      ax.set_ylabel(r"$y / D_{j}$", fontsize = 16)
      ax.set_yticks(np.array([0, 5, 10]))
    else:
      ax.set_yticks([])
plt.savefig("zone_division.png", dpi=300, bbox_inches="tight")
#%%
ncols = 2; nrows = 2
field_name = "zone"
vmin=1.0; vmax = 4.0
cmap = "jet"
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(4, 6))
k0 = 1.5
xl = np.linspace(-4, 0, 50)
xr = np.linspace( 0, 4, 50)
for ipx in range(0, ncols):
  for ipy in range(0, nrows):
    ix = ipx * ncols + ipy
    ax = axs[ipx,ipy]
    phi = tet_xs[ix][field_name]
    mf = tet_xs[ix]["mixture_fraction"]
    im = ax.imshow(phi, origin="lower", aspect='equal',
                vmin=vmin, vmax=vmax, 
                cmap=cmap, extent=extent_x)
    ax.contour(mf, levels=[0.0252], aspect='equal',
          origin='lower', 
          colors=['black'], extent=extent_x)
    xloc = xps[ix]
    zloc = get_trajectory_z(xloc)

    if (False):
      y0 = -2 
      z0 = zloc / Djet
      b0 = z0 - k0 * y0 
      ax.plot(xl, k0*xl+b0,color="cyan", linestyle="--")
      y0 = -2 
      z0 = zloc / Djet
      b0 = z0 - (-k0) * y0 
      ax.plot(xl, (-k0)*xl+b0,color="cyan", linestyle="--")
      y0 = 2 
      z0 = zloc / Djet
      b0 = z0 - k0 * y0 
      ax.plot(xr, k0*xr+b0,color="cyan", linestyle="--")
      y0 = 2 
      z0 = zloc / Djet
      b0 = z0 - (-k0) * y0 
      ax.plot(xr, (-k0)*xr+b0,color="cyan", linestyle="--")
      ax.set_ylim([0, 12])
      ax.scatter(-2, zloc/Djet, color="magenta")
      ax.scatter(+2, zloc/Djet, color="magenta")
    ax.set_title(r"$x/D_j=" + str(int(xloc/Djet))+"$")

    if ipx == 0 and ipy == 0:
      cax = ax.inset_axes([-3, 10.0, 6, 1.0], transform=ax.transData)
      cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                        ticks=[1.5, 2.5, 3.5])
      cb.ax.xaxis.set_tick_params(color="black")
      cb.ax.yaxis.set_tick_params(color="black")
      cb.outline.set_edgecolor("black")
      cb.ax.set_xticklabels(["Lee", "Wind", "Int."], color="black", fontsize=10)  # horizontal colorbar

    if False:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes('right', size='5%', pad=0.05)
      fig.colorbar(im, cax=cax, orientation='vertical')
    if (ipx == nrows-1):
      ax.set_xlabel(r"$z / D_\mathrm{j}$", fontsize = 14)
      ax.set_xticks(np.array([-4, -2, 0, 2, 4]))
    else:
        #ax.set_xlabel([])
      ax.set_xticks(np.array([]))

    if (ipy == 0):
      ax.set_ylabel(r"$y / D_{j}$", fontsize = 14)
      ax.set_yticks(np.array([0, 5, 10]))
    else:
      ax.set_yticks([])
   

#%%