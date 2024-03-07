#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
import h5py

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive/plt_08800"
ds = yt.load(fn)
dd = ds.all_data()
zst = 0.0252
#%%
xmin = -17.5E-4; xmax = 122.5E-4
ymin = -2E-3; ymax = 2E-3
zmin = -0.0E-3; zmax = 6.0E-3
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

# Normalized parameters
tau_flow = Ly / 20. 

res_x=1000
res_z=res_x * (Lz / Lx)
res_y=res_x * (Ly / Lx)
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels
res_yz = [res_y, res_z]  # create an image with 1000x1000 pixels


loc = 2.0*5E-4
loc = 2.2188e+00*5E-4
ye = 1.5625e-01 * 5E-4
ze = 1.5938e+00 * 5E-4

fn = "temp"; dir = "y"; vmin=400; vmax=750
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
arr_mf = np.array(frb["mixture_fraction"])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
ax.scatter(ye, ze, s=15, color="r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$T \; [\mathrm{K}]$", fontsize = 20)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")
#%%
fn = "mixture_fraction"; dir = "y"; vmin = 0.0; vmax = 0.5
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax])
ctr = ax.contour(arr_temp.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "mixture_fraction"; dir = "y"; vmin = 0.0; vmax = 1.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax])
ctr = ax.contour(arr_temp.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
#plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

field_name = "mixture_fraction"
fig, ax = plt.subplots()
spl = yt.SlicePlot(ds, "y", 
             [(field_name)], 
             center=(xmin+0.5*Lx, loc,zmin+0.5*Lz), width=((Lz, "cm"),(Lx, "cm")),)
spl.set_log((field_name), False)
spl.annotate_grids()
#%%
fn = "HeatRelease"; dir = "y"; vmin = 0; vmax = 1E11
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
phi = np.log10(np.abs(arr_temp.transpose()))
phi = arr_temp.transpose()
im = ax.imshow(phi, vmin=vmin, vmax=vmax, origin="lower", cmap="hot", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "x_velocity"; dir = "y"; vmin = -120; vmax = 120
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="seismic", extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "Y(H2)"; dir = "x"; loc = 1.0E-3
slc = ds.slice(dir, coord=loc, center=(loc, ymin+0.5*Ly, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Ly, "cm"),(Lz, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[zmin, zmax, ymin, ymax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "temp"; dir = "x"; loc = xmax-0.001
slc = ds.slice(dir, coord=loc, center=(loc, ymin+0.5*Ly, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Ly, "cm"),(Lz, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[zmin, zmax, ymin, ymax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "mixture_fraction"; dir = "z"; loc = 0.002
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[ymin, ymax, xmin, xmax])
ctr = ax.contour(arr_temp.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[ymin, ymax, xmin, xmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "density"; dir = "z"; loc = 0.00001
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[ymin, ymax, xmin, xmax])
ctr = ax.contour(arr_temp.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[ymin, ymax, xmin, xmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")
#%%
