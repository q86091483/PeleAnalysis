#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
import yt
import h5py
from copy import copy

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_reactive/plt_05000"
fn = "/scratch/b/bsavard/zisen347/scopingRuns/MicroMix/plt_00000"

ds = yt.load(fn) # unit_system="cgs")
zst = 0.0252
#%%
xmin = -17.5E-4; xmax = 122.5E-4
#xmin = -7.5E-4; xmax = 1.325E-2
ymin = -2E-3; ymax = 2E-3
zmin = -0.0E-3; zmax = 6.0E-3
Lx = xmax - xmin; nx0 = 224
Ly = ymax - ymin; ny0 = 128
Lz = zmax - zmin; nz0 = 96

# Normalized parameters
tau_flow = Ly / 20. 

finest_level = 3;
res_x=nx0 * 2**finest_level
res_z=ny0 * 2**finest_level
res_y=nz0 * 2**finest_level
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels
res_yz = [res_y, res_z]  # create an image with 1000x1000 pixels

coord = None #(-1.5625e-01,-1.9688e+00,1.2188e+00)
lref = 5.0E-4
if coord != None:
  xe = coord[0] * lref
  loc = coord[1]  * lref
  ze = coord[2] * lref
else:
  xe = 0.0 * lref
  loc = 0.0 * lref
  ze = 0.0 * lref

fn = "temp"; dir = "y"; vmin = 600; vmax = 900
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")


fn = "x_velocity"; dir = "y"; vmin = 0.0; vmax = 25
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
ax.set_title(r"$ U $", fontsize=18)
ax.set_xlabel(r"$x$", fontsize=18)
ax.set_xlabel(r"$z$", fontsize=18)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+"_interp.png", dpi=300, bbox_inches="tight")
#%%
fn = "y_velocity"; dir = "y"; vmin = -4.0; vmax = 4.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
ax.set_title(r"$ V $", fontsize=18)
ax.set_xlabel(r"$x$", fontsize=18)
ax.set_xlabel(r"$z$", fontsize=18)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+"_interp.png", dpi=300, bbox_inches="tight")
#%%
fn = "z_velocity"; dir = "y"; vmin = -4.0; vmax = 4.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
ax.set_title(r"$ W $", fontsize=18)
ax.set_xlabel(r"$x$", fontsize=18)
ax.set_xlabel(r"$z$", fontsize=18)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+"_interp.png", dpi=300, bbox_inches="tight")


#%%