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

case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/"
fns = glob.glob(case_folder + "*.h5")

xmin = -17.5E-4; xmax = 122.5E-4
ymin = -2.0E-3; ymax = 2.0E-3
zmin = -0.0E-3; zmax = 6.0E-3
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin
#fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt_11000_Y.h5"
fn = fns[0]
print(fn)
f = h5py.File(fn, 'r+')
wt_sum = np.array(f["DATA"]["volume_mean"])
rho_wtsum = np.array(f["DATA"]["rho_mean"]) 
ts11_wtsum = np.array(f["DATA"]["ts11_mean"]) 
ts22_wtsum = np.array(f["DATA"]["ts22_mean"]) 
ts33_wtsum = np.array(f["DATA"]["ts33_mean"]) 
ts12_wtsum = np.array(f["DATA"]["ts12_mean"]) 
ts13_wtsum = np.array(f["DATA"]["ts13_mean"]) 
ts23_wtsum = np.array(f["DATA"]["ts23_mean"]) 
mu_wtsum = np.array(f["DATA"]["mu_mean"]) 
ts_wtsum = np.array(f["DATA"]["ts_sum_mean"])

rho = rho_wtsum / wt_sum
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
fig, ax = plt.subplots()
vmin = 1E-7; vmax = 1E-2
im = ax.imshow(eta[:,5,:].transpose(), cmap="tab20c", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
#%%
for ifn, fn in enumerate(fns[1:]):
  f = h5py.File(fn, 'r+')
  wt_sum = wt_sum + np.array(f["DATA"]["volume_mean"])
  rho_wtsum = rho_wtsum + np.array(f["DATA"]["rho_mean"]) 
  ts11_wtsum = ts11_wtsum + np.array(f["DATA"]["ts11_mean"]) 
  ts22_wtsum = ts22_wtsum + np.array(f["DATA"]["ts22_mean"]) 
  ts33_wtsum = ts33_wtsum + np.array(f["DATA"]["ts33_mean"]) 
  ts12_wtsum = ts12_wtsum + np.array(f["DATA"]["ts12_mean"]) 
  ts13_wtsum = ts13_wtsum + np.array(f["DATA"]["ts13_mean"]) 
  ts23_wtsum = ts23_wtsum + np.array(f["DATA"]["ts23_mean"]) 
  mu_wtsum = mu_wtsum + np.array(f["DATA"]["mu_mean"]) 
  ts_wtsum = ts_wtsum + np.array(f["DATA"]["ts_sum_mean"]) 

rho = rho_wtsum / wt_sum
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
fig, ax = plt.subplots()
vmin = 1E-7; vmax = 1E-2
im = ax.imshow(eta[:,5,:].transpose(), cmap="tab20c", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig("eta.png", dpi=300, bbox_inches="tight")
#%%