#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

import cantera as ct
gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/nuig_H2_4atm/chem.yaml")
Nsp = len(gas_mix.species_names)
H2_ID = gas_mix.species_index("H2")
O2_ID = gas_mix.species_index("O2")
N2_ID = gas_mix.species_index("N2")
NO_ID = gas_mix.species_index("NO")


data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/"
fn_derived = os.path.join(data_folder, "Slice2D_derived/Micromix/Derived_y=-1.125E-03/plt_10580_derived_t=6.077E-04.npz")
fn_plt     = os.path.join(data_folder, "Slice2D_plt_lev1/Micromix/HRR_T_y=-1.125E-03/plt_10580_t=6.077E-04.npz")
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3
Djet = 4.5E-04
extent=np.array([xmin, xmax, zmin, zmax]) / Djet

# Read data
P0      = 101325*4.0
tet_plt = np.load(fn_plt)
tet_der = np.load(fn_derived)
temp    = tet_plt["temp"]
density = tet_plt["density"]
pres    = np.zeros_like(density) + P0
hrr     = tet_plt["HeatRelease"]
hrrFI   = tet_der["HeatReleaseFI"]
mixfrac = tet_der["mixture_fraction"]
pv      = tet_der["pv"]
rho     = tet_der["rho"]
FI      = tet_der["FI"]
print("Max/min rho: ",  np.amin(rho), np.amax(rho))
print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min progvar: ",  np.amin(pv), np.amax(pv))
print("Max/min FI: ",  np.amin(FI), np.amax(FI))
nx = temp.shape[0] 
ny = temp.shape[1]

tet_Y = np.zeros((nx, ny, Nsp))
for isp in range(0, Nsp):
  spn = gas_mix.species_names[isp]
  tet_Y[:,:,isp] = tet_plt["Y("+spn+")"][:,:]
#%%
fstate = ct.SolutionArray(gas_mix, (nx, ny))

fstate[0].TPY = temp[0,:], pres[0,:], tet_Y[0,:,:]

#%%
for i in range(0, nx):
  fstate[i].TPY = temp[i,:], pres[i,:], tet_Y[i,:,:]



#%%
# T in state
fig, ax = plt.subplots()
vmin = 300; vmax = 2500
im = ax.imshow(fstate.T[:,:], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="jet", extent=extent)
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
vmin = 300; vmax = 2500
im = ax.imshow(temp, origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="jet", extent=extent)
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
#%%

# Y(NO) in state
fig, ax = plt.subplots()
vmin = 0; vmax = 1E-4
im = ax.imshow(fstate.Y[:,:,NO_ID], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="jet", extent=extent)
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
vmin = 0; vmax = 1E-4
im = ax.imshow(tet_Y[:,:,NO_ID], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="jet", extent=extent)
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
#%%