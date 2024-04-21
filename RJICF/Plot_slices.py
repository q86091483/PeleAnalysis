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
#fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive_1/plt_04370"
#fn = "/scratch/b/bsavard/zisen347/scopingRuns/NUIG_Re3600_2J6_4atm/plt_02300"
fn = "/scratch/b/bsavard/zisen347/scopingRuns/MicroMix/plt_09620"

ds = yt.load(fn) # unit_system="cgs")
zst = 0.0252
#%%
xmin = -15.75E-4; xmax = 112.25E-4
#xmin = -7.5E-4; xmax = 1.325E-2
ymin = -1.8E-3; ymax = 1.8E-3
zmin = -0.0E-3; zmax = 5.6E-3
Lx = xmax - xmin; nx0 = 256
Ly = ymax - ymin; ny0 = 72
Lz = zmax - zmin; nz0 = 112

# Normalized parameters
tau_flow = Ly / 20. 

finest_level = 2;
res_x=nx0 * 2**finest_level
res_z=ny0 * 2**finest_level
res_y=nz0 * 2**finest_level
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels
res_yz = [res_y, res_z]  # create an image with 1000x1000 pixels

coord = None #(-1.5625e-01,-1.9688e+00,1.2188e+00)
lref = 0.45E-3
if coord != None:
  xe = coord[0] * lref
  loc = coord[1]  * lref
  ze = coord[2] * lref
else:
  xe = 2.0 * lref
  loc = -2.0 * lref
  ze = 2.0 * lref

fn = "mixture_fraction"; dir = "y"; vmin=0.0; vmax=1.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_mf = np.array(frb["mixture_fraction"])
fig, ax = plt.subplots()
im = ax.imshow(arr_mf.transpose(), 
              origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
ax.scatter(xe, ze, s=15, color="r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$T \; [\mathrm{K}]$", fontsize = 18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "HeatRelease"; dir = "y"; vmin = 0; vmax = 1E11
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
phi = np.log10(np.abs(arr_temp.transpose()))
phi = arr_temp.transpose()
im = ax.imshow(phi, vmin=vmin, vmax=vmax, origin="lower", cmap="hot", extent=[xmin, xmax, zmin, zmax])
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$\mathrm{HRR}$", fontsize = 18)
ax.set_xlabel(r"$x \; [\mathrm{m}$]", fontsize=18)
ax.set_ylabel(r"$y \; [\mathrm{m}$]", fontsize=18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "temp"; dir = "y"; vmin=400; vmax=2600
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
arr_mf = np.array(frb["mixture_fraction"])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
ax.scatter(xe, ze, s=15, color="r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$T \; [\mathrm{K}]$", fontsize = 18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "Y(OH)"; dir = "y"; vmin=0.0; vmax=2E-2
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
ax.scatter(xe, ze, s=15, color="r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$Y_{OH}$", fontsize = 18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "density"; dir = "y"; vmin = 0.0; vmax = 5.0
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax])
ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title(r"$\rho \; [\mathrm{kg/m^3}]$", fontsize = 18)
ax.set_xlabel(r"$x \; [\mathrm{m}]$", fontsize=18)
ax.set_ylabel(r"$y \; [\mathrm{m}]$", fontsize=18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "Y(NO)"; dir = "y"; vmin = 0.0; vmax = 1E-4
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax])
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title(r"$Y_\mathrm{NO}$", fontsize = 18)
ax.set_xlabel(r"$x \; [\mathrm{m}]$", fontsize=18)
ax.set_ylabel(r"$y \; [\mathrm{m}]$", fontsize=18)
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")
#%%



#%%
#fig = plt.figure()
fn  = "mag_vort"
spl = yt.SlicePlot(ds, "y", 
             [(fn)], 
             origin='native',
             center=(xmin+0.5*Lx, loc,zmin+0.5*Lz), width=((Lz, "cm"),(Lx, "cm")), )
spl.set_log((fn), False)
spl.annotate_grids(cmap='brg_r')
spl.set_cmap('mag_vort', 'binary')
spl.set_zlim('mag_vort', 0, 1e6)
spl.set_xlabel(r"$x \; \mathrm{[m]}$")
spl.set_ylabel(r"$z \; \mathrm{[m]}$")
spl.swap_axes()
spl.save()

#%%
fn  = "HeatRelease"
spl = yt.SlicePlot(ds, "y", 
             [(fn)], 
             origin='native',
             center=(xmin+0.5*Lx, loc,zmin+0.5*Lz), width=((Lz, "cm"),(Lx, "cm")), )
spl.set_log((fn), False)
spl.annotate_grids(cmap='brg_r')
spl.set_cmap(fn, 'binary')
spl.set_zlim(fn, 0, 1e11)
spl.set_xlabel(r"$x \; \mathrm{[m]}$")
spl.set_ylabel(r"$z \; \mathrm{[m]}$")
spl.swap_axes()
spl.save()

#%%

fn = "mag_vort"; dir = "y"; vmin = 0; vmax = 1E6
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
phi = np.log10(np.abs(arr_temp.transpose()))
phi = arr_temp.transpose()
im = ax.imshow(phi, vmin=vmin, vmax=vmax, origin="lower", cmap="binary", extent=[xmin, xmax, zmin, zmax])
ctr = ax.contour(arr_mf, levels=[zst] ,origin='lower', colors=['white'], extent=[xmin, xmax, zmin, zmax],
                 linestyle = "--")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%
fn = "x_velocity"; dir = "y"; vmin = -50; vmax = 200
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), vmin=vmin, vmax=vmax, origin="lower", cmap="jet", extent=[xmin, xmax, zmin, zmax])
ax.scatter(xe, ze, s=15, color="g")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%
height = 0.002
fn = "temp"; dir = "z"; loc = height
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[ymin, ymax, xmin, xmax])
arr_mf = np.array(frb[fn])
ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[ymin, ymax, xmin, xmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

fn = "HeatRelease"; dir = "z"; loc = height
vmin = 0; vmax=2E11
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="hot", 
               vmin=vmin, vmax=vmax,
               extent=[ymin, ymax, xmin, xmax])
arr_mf = np.array(frb[fn])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")


#%%
height = 0*lref
fn = "Y(H2)"; dir = "z"; loc = height
vmin = 0.0; vmax = 1.0
palette = copy(plt.get_cmap('jet_r'))
palette.set_over('white', 0.0) 
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, ymin+0.5*Ly, loc), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res_xy)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), origin="lower", cmap="viridis", extent=[ymin, ymax, xmin, xmax],
               vmin = vmin, vmax = vmax)
#ax.add_patch(plt.Circle((1*lref, 0.0), lref, color='b', fill=False, linestyle="--"))              
#ax.add_patch(plt.Circle((1*lref, 0.0), lref, color='b', fill=False, linestyle="--"))              
arr_mf = np.array(frb[fn])
#ctr = ax.contour(arr_mf.transpose(), levels=[zst] ,origin='lower', colors=['white'], extent=[ymin, ymax, xmin, xmax])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")

#%%