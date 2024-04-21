#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import glob as glob
import os

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Input
# Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_FlameStructure_basic/"
# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_070*_Y.h5")
# Domain length of Pele Case
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = 1.8E-3
zmin = -0.0E-3; zmax = 5.6E-3
nx = 128; 
ny = 36;
nz = 56;
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

#
A = 0.4
B = 0.28
J = 6.66
Djet = 4.5E-4
xs = np.linspace(0, 0.01, 1000)
ys = J * Djet * A * np.power(xs / (J*Djet), B)
#%%
# Declaration
wt_sum = np.zeros((nx, ny, nz))
rho_wtsum = np.zeros((nx, ny, nz))
T_wtsum = np.zeros((nx, ny, nz))
HRR_wtsum = np.zeros((nx, ny, nz))
ts11_wtsum = np.zeros((nx, ny, nz))
ts22_wtsum = np.zeros((nx, ny, nz))
ts33_wtsum = np.zeros((nx, ny, nz))
ts12_wtsum = np.zeros((nx, ny, nz))
ts13_wtsum = np.zeros((nx, ny, nz))
ts23_wtsum = np.zeros((nx, ny, nz))
mu_wtsum = np.zeros((nx, ny, nz))
ts_wtsum = np.zeros((nx, ny, nz))

for ifn, fn in enumerate(fns[0:]):
  f = h5py.File(fn, 'r+')
  wt_sum = wt_sum + np.array(f["DATA"]["volume_mean"])
  rho_wtsum = rho_wtsum + np.array(f["DATA"]["rho_mean"]) 
  T_wtsum = T_wtsum + np.array(f["DATA"]["temp_mean"]) 
  HRR_wtsum = HRR_wtsum + np.array(f["DATA"]["HeatRelease_mean"]) 
  ts11_wtsum = ts11_wtsum + np.array(f["DATA"]["ts11_mean"]) 
  ts22_wtsum = ts22_wtsum + np.array(f["DATA"]["ts22_mean"]) 
  ts33_wtsum = ts33_wtsum + np.array(f["DATA"]["ts33_mean"]) 
  ts12_wtsum = ts12_wtsum + np.array(f["DATA"]["ts12_mean"]) 
  ts13_wtsum = ts13_wtsum + np.array(f["DATA"]["ts13_mean"]) 
  ts23_wtsum = ts23_wtsum + np.array(f["DATA"]["ts23_mean"]) 
  mu_wtsum = mu_wtsum + np.array(f["DATA"]["mu_mean"]) 
  ts_wtsum = ts_wtsum + np.array(f["DATA"]["ts_sum_mean"]) 

rho = rho_wtsum / wt_sum
T = T_wtsum / wt_sum
HRR = HRR_wtsum / wt_sum
ts11 = ts11_wtsum / wt_sum
ts22 = ts22_wtsum / wt_sum
ts33 = ts33_wtsum / wt_sum
ts12 = ts12_wtsum / wt_sum
ts13 = ts13_wtsum / wt_sum
ts23 = ts23_wtsum / wt_sum
mu = mu_wtsum / wt_sum
ts = ts_wtsum / wt_sum
sum_of_meanprod = ts11*ts11 + ts22*ts22 + ts33*ts33+ \
                 2*ts12*ts12 + 2*ts13*ts13 + 2*ts23*ts23
eps = ts - sum_of_meanprod
eps = 2 * eps * mu / rho
eta = np.power(mu/rho, 3) / eps
eta = np.power(eta, 0.25)

# T - y
fig, ax = plt.subplots()
vmin = 500; vmax = 2800
im = ax.imshow(T[:,32,:].transpose(), cmap="jet", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
ax.plot(xs, ys)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$ T \; [\mathrm{K}]$", fontsize=20)
plt.savefig("T.png", dpi=300, bbox_inches="tight")
#%%
# HRR - y
fig, ax = plt.subplots()
vmin = 0; vmax = 1E11
phi = HRR[:,27,:]
im = ax.imshow(phi.transpose(), cmap="hot", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$HRR \; [\mathrm{K}]$", fontsize=20)
plt.savefig("T.png", dpi=300, bbox_inches="tight")

#%%
fig, ax = plt.subplots()
vmin = 0; vmax = 5
im = ax.imshow(rho[:,8,:].transpose(), cmap="jet", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig("eta.png", dpi=300, bbox_inches="tight")#
#%%