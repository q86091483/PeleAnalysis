#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import yt

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt_11000_derived"
#fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive/plt_11000"
ds = yt.load(fn)
dd = ds.all_data()
#%%
xmin = -17.5E-4; xmax = 122.5E-4
ymin = -2.0E-3; ymax = 2.0E-3
zmin = -0.0E-3; zmax = 6.0E-3
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

res_x=500
res_z=res_x * (Lz / Lx)
res_y=res_x * (Ly / Lx)
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels

loc = 5E-4 * 2
fn = "rho"; dir = "y"; vmin=1E-7; vmax=2E-4
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
#arr_temp = np.log10(np.abs(np.array(frb[fn])))
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="tab20", extent=[xmin, xmax, zmin, zmax],
              norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%
fn = "rho"; dir = "y"; vmin=0.0; vmax=5.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=0.0, vmax=vmax,
              origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax],
              )
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%
fn = "ts23"; dir = "y"; 
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), # vmin=0.0, vmax=vmax,
              origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax],
              )
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%