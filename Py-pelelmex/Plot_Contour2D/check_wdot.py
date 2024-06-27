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
#gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/nuig_H2_4atm/chem.yaml")
gas_mix = ct.Solution("/home/b/bsavard/zisen347/codes/exascale/PelePhysics/Support/Mechanism/Models/nuig_H2_4atm_CHNO/nuig_H2_32sp_4atm.yaml")
Nsp = len(gas_mix.species_names)
H2_ID = gas_mix.species_index("H2")
O2_ID = gas_mix.species_index("O2")
N2_ID = gas_mix.species_index("N2")
NO_ID = gas_mix.species_index("NO")
N2O_ID = gas_mix.species_index("N2O")
NNH_ID = gas_mix.species_index("NNH")

data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/"
fn_derived = os.path.join(data_folder, "Slice2D_derived/Micromix/Derived_y=-1.125E-03/plt_10580_derived_t=6.077E-04.npz")
fn_plt     = os.path.join(data_folder, "Slice2D_plt_lev1/Micromix/HRR_T_y=-1.125E-03/plt_10580_t=6.077E-04.npz")
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3
Djet = 4.5E-04
extent=np.array([xmin, xmax, zmin, zmax]) / Djet
#%%
# Read data
P0      = 101325*4.0
tet_plt = np.load(fn_plt)
tet_der = np.load(fn_derived)
temp    = tet_der["temp"]
density = tet_plt["density"]
pres    = np.zeros_like(density) + P0
hrr     = tet_plt["HeatRelease"]
hrrFI   = tet_der["HeatReleaseFI"]
mixfrac = tet_der["mixture_fraction"]
pv      = tet_der["pv"]
rho     = tet_der["rho"]
rhorr_NO = tet_der["rhorr(NO)"]
rhorr_N2O = tet_der["rhorr(N2O)"]
rhorr_NNH = tet_der["rhorr(NNH)"]


FI      = tet_der["FI"]
R10     = tet_der["R10"]
print("Max/min rho: ",  np.amin(rho), np.amax(rho))
print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min progvar: ",  np.amin(pv), np.amax(pv))
print("Max/min FI: ",  np.amin(FI), np.amax(FI))
nx = temp.shape[0] 
ny = temp.shape[1]

# Array of species
tet_Y = np.zeros((nx, ny, Nsp))
for isp in range(0, Nsp):
  spn = gas_mix.species_names[isp]
  tet_Y[:,:,isp] = tet_plt["Y("+spn+")"][:,:]

# Load data into ct.SolutionArray
fstate = ct.SolutionArray(gas_mix, (nx, ny))
for i in range(0, nx):
  fstate[i].TPY = temp[i,:], pres[i,:], tet_Y[i,:,:]

# Way 1
net_rates_of_progress = fstate.net_rates_of_progress
forward_rates_of_progress = fstate.forward_rates_of_progress
reverse_rates_of_progress = fstate.reverse_rates_of_progress

net_rates_of_progress_ddC = fstate.net_rates_of_progress_ddC
delta_enthalpy = fstate.delta_enthalpy
hrr_state = np.sum(net_rates_of_progress * delta_enthalpy, axis=(2))
# Way 2
net_production_rates = fstate.net_production_rates
partial_molar_enthalpies = fstate.partial_molar_enthalpies
hrr_state = np.sum(net_production_rates * partial_molar_enthalpies, axis=(2))

#%%
# Reaction fstate
fig, ax = plt.subplots()
vmin = -10; vmax = 10
im = ax.imshow(R10[:,:]/1000, origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="seismic", extent=extent)
ax.set_title(r"$Pele: " + str(gas_mix.reactions()[22]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], ext0ent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
vmin = -10; vmax = 10
im = ax.imshow(fstate.net_rates_of_progress[:,:,22], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="seismic", extent=extent)
ax.set_title(r"$Cantera: " + str(gas_mix.reactions()[22]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
vmin = -100; vmax = 100
im = ax.imshow(R10[:,:]/1000-fstate.net_rates_of_progress[:,:,22], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="seismic", extent=extent)
ax.set_title(r"$Pele-Cantera: " + str(gas_mix.reactions()[22]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower', 
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')


#%%
if (True):
  # rhorr(NO) in fstate
  fig, ax = plt.subplots()
  vmin = -0.1; vmax = 0.1
  im = ax.imshow(net_production_rates[:,:,N2O_ID]*gas_mix.molecular_weights[NO_ID], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="seismic", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')

  fig, ax = plt.subplots()
  vmin = -0.1; vmax = 0.1
  im = ax.imshow(rhorr_N2O, origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="seismic", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')


  fig, ax = plt.subplots()
  vmin = 0; vmax = 1E11
  im = ax.imshow(-hrr_state[:,:], origin="lower", 
               vmin=vmin, vmax=vmax, 
                cmap="hot", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')

#%%