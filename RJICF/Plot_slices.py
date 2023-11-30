#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
fn = "/home/zisen347/scratch/scoping_runs/exascale/MicroMix/plt_01200"
ds = yt.load(fn)
dd = ds.all_data()

xmin = -7.5E-4; xmax = 11.25E-3;
ymin = -1.5E-3; ymax = 1.5E-3
zmin = 0; zmax = 2.8E-3
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

res_x=500
res_z=res_x * (Lz / Lx)
res_y=res_x * (Ly / Lx)
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels

loc = 7.5E-4
fn = "temp"; dir = "y"; vmin=400; vmax=2400
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "Y(H2)"; dir = "y"; vmin = 0.0; vmax = 1.E0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "Y(H2)"; dir = "x"; loc = 7.5E-4
slc = ds.slice(dir, coord=loc, center=(loc, ymin+0.5*Ly, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Ly, "cm"),(Lz, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[zmin, zmax, ymin, ymax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "temp"; dir = "z"; loc = 0.0000
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[ymin, ymax, xmin, xmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%