#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import glob as glob
import os
import os.path as path
import re
import pandas as pd
import cantera as ct

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Input
# Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix/"
# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_07*_Y.h5")
# Data folder - where 1D flame data is stored
data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/UnstrainedPremixed1D"

zout = np.linspace(0, 1, 51)
pvout = np.linspace(0, 1, 51)

# Domain length of Pele Case
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = 1.8E-3
zmin = -0.0E-3; zmax = 5.6E-3
nmix = 30;
#%%
# Files to read
# Which flame.csv to process
csv_pattern = "Z=*.csv"
fns1D_unsorted = glob.glob(path.join(data_folder, csv_pattern))
def get_Z_key(s):
  ss = re.split("/", s)[-1]
  ss = re.split("=", ss)[-1]
  ss = re.split(".csv", ss)[0]
  return(float(ss))
fns1D_sorted = sorted(fns1D_unsorted, key=get_Z_key)
Zs1D = []
for ifn, fn in enumerate(fns1D_sorted):
  Zs1D.append(get_Z_key(fn))
#%%
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

# Read parameters
fn = fns[0]
f = h5py.File(fn, 'r+')
tet = np.array(f["DATA"]["volume_mean"])
nx = tet.shape[0]
ny = tet.shape[1]
nz = tet.shape[2]
nmf = tet.shape[3]
npv = tet.shape[4]

#%%
mech = "nuig_H2_4atm/chem.yaml"
def get_states(zs, equilibrate = True):
  gas_f = ct.Solution(mech)
  gas_o = ct.Solution(mech)
  gas_m = ct.Solution(mech)
  states = ct.SolutionArray(gas_m) 
  for iz, z in enumerate(zs):
    # Fuel and oxidizer stream
    X_f    = {}; X_f["H2"] = 1.0; X_f["N2"] = 1 - X_f["H2"] 
    X_o    = {}; X_o["O2"] = 0.21; X_o["N2"] = 0.79
    gas_f.TPX = 300, 405300, X_f
    gas_o.TPX = 750, 405300, X_o
    # Mixture with z
    Ym = z*gas_f.Y + (1-z)*gas_o.Y
    Hm = z*gas_f.enthalpy_mass + (1-z)*gas_o.enthalpy_mass
    Pm = gas_f.P
    gas_m.HPY = Hm, Pm, Ym
    states.append(T = gas_m.T,
                  P = gas_m.P,
                  Y = gas_m.Y)
    if equilibrate == True:
      states.equilibrate("HP")
  return states

states_0D_min = get_states(zout, equilibrate=False)
states_0D_max = get_states(zout, equilibrate=True)

#%%
wt_sum = np.zeros((nx, ny, nz, nmf, npv))
T_wtsum = np.zeros((nx, ny, nz, nmf, npv))
mf_wtsum = np.zeros((nx, ny, nz, nmf, npv))
pv_wtsum = np.zeros((nx, ny, nz, nmf, npv))
hrr_wtsum = np.zeros((nx, ny, nz, nmf, npv))
rhorr_NO_wtsum = np.zeros((nx, ny, nz, nmf, npv))
rhorr_N2O_wtsum = np.zeros((nx, ny, nz, nmf, npv))
rhorr_NNH_wtsum = np.zeros((nx, ny, nz, nmf, npv))

wt_xyzZC_sum = wt_sum + f["DATA"]["volume_mean"]
rho_xyzZC_sum = wt_sum + f["DATA"]["rho_mean"]
T_xyzZC_wtsum = T_wtsum + f["DATA"]["temp_mean"]
mf_xyzZC_wtsum = mf_wtsum + f["DATA"]["mixture_fraction_mean"]
pv_xyzZC_wtsum = pv_wtsum + f["DATA"]["pv_mean"]
hrr_xyzZC_wtsum = hrr_wtsum + f["DATA"]["HeatRelease_mean"]
rhorr_NNH_xyzZC_wtsum = rhorr_NNH_wtsum + f["DATA"]["rhorr(NNH)_mean"]
rhorr_N2O_xyzZC_wtsum = rhorr_N2O_wtsum + f["DATA"]["rhorr(N2O)_mean"]
rhorr_NO_xyzZC_wtsum = rhorr_NO_wtsum + f["DATA"]["rhorr(NO)_mean"]

#%%

axis_sum = (0, 1, 2)
wt_ZC_sum = np.sum(wt_xyzZC_sum, axis=axis_sum)
T_ZC_sum = np.sum(T_xyzZC_wtsum, axis=axis_sum)
mf_ZC_sum = np.sum(mf_xyzZC_wtsum, axis=axis_sum)
pv_ZC_sum = np.sum(pv_xyzZC_wtsum, axis=axis_sum)
hrr_ZC_sum = np.sum(hrr_xyzZC_wtsum, axis=axis_sum)
rhorr_NO_ZC_sum = np.sum(rhorr_NO_xyzZC_wtsum, axis=axis_sum)
rhorr_N2O_ZC_sum = np.sum(rhorr_N2O_xyzZC_wtsum, axis=axis_sum)
rhorr_NNH_ZC_sum = np.sum(rhorr_NNH_xyzZC_wtsum, axis=axis_sum)

mf_ZC = mf_ZC_sum / wt_ZC_sum
pv_ZC = pv_ZC_sum / wt_ZC_sum
T_ZC = T_ZC_sum / wt_ZC_sum
hrr_ZC = hrr_ZC_sum / wt_ZC_sum
rhorr_NO_ZC = rhorr_NO_ZC_sum / wt_ZC_sum
rhorr_N2O_ZC = rhorr_N2O_ZC_sum / wt_ZC_sum
rhorr_NNH_ZC = rhorr_NNH_ZC_sum / wt_ZC_sum

figsize = (5,4)
labelsize = 20
#%% T
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(T_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          vmin=300, vmax=2500, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle T \; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)

#%% HRR
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(hrr_ZC)
phi = np.log10(phi)
ax.imshow(phi.transpose(), origin="lower", 
          vmin=6, vmax=12, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle \mathrm{log}_{10}(\mathrm{HRR}) \; | \; Z, C \rangle$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)

#%% C
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(pv_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          vmin=0, vmax=1, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle C \; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)

#%% Z
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(mf_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          vmin=0, vmax=1, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle Z \; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)

#%% rhorr_NO
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(rhorr_NO_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          # vmin=0, vmax=5.0, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle \rho \dot{\omega}_\mathrm{NO}\; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)
#%% rhorr_NNH
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(rhorr_NNH_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          # vmin=0, vmax=5.0, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle \rho \dot{\omega}_\mathrm{NNH}\; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)
#%% rhorr_N2O
fig, ax = plt.subplots(figsize=figsize)
phi = np.abs(rhorr_N2O_ZC)
ax.imshow(phi.transpose(), origin="lower", 
          # vmin=0, vmax=5.0, 
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.0])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle \rho \dot{\omega}_\mathrm{N2O}\; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize-2)
ax.tick_params(axis='both', which='minor', labelsize=labelsize-2)




# %%
