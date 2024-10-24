#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
#import matplotlib as mpl
#mpl.rcParams.update(mpl.rcParamsDefault)
import glob as glob
import os
import os.path as path
import re
import pandas as pd
import cantera as ct
from mixture_fraction import mf

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#%%

# Input
# Where to output the result
#case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/backup_res_MicroMix_xzone"
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix_xzone_lev3"

# Where to read h5 files that contain condition smean
fns = glob.glob(os.path.join(case_folder, "plt_*.h5"))
# Data folder - where 1D flame data is stored
data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/UnstrainedPremixed1D"
# Field names
field_names = [["T", "wdot(O2)", "wdot(H2)", "wdot(H2O)"],
                ["wdot(NH)", "wdot(NNH)", "wdot(OH)", "wdot(H2O2)"],
               ["wdot(NO)", "wdot(N2O)", "wdot(NO2)", "wdot(N)"],
              ]
# Derived parameter
fn = fns[0]
f = h5py.File(fn, 'r+')
tet = np.array(f["DATA"]["volume_mean"])
nx = tet.shape[0]
ny = tet.shape[1]
nz = tet.shape[2]
nmf = tet.shape[3]
npv = tet.shape[4]
nzo = tet.shape[5]
nfi = tet.shape[6]
zout = np.linspace(0, 1, nmf)
pvout = np.linspace(0, 1, npv)
labelsize = 16
# Domain length of Pele Case
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = 1.8E-3
zmin = -0.0E-3; zmax = 5.6E-3
Lx = xmax - xmin
Ly = ymax - ymin
Lz = zmax - zmin
nmix = 30;
# Mechanism
mech = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/nuig_H2_4atm/chem.yaml"
gas1D = ct.Solution(mech)
Nsp = len(gas1D.species_names)
species_names = gas1D.species_names
#%%
# Load 1D premixed flames
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
Zs1D = np.array(Zs1D)

# Load 0D equilibirum state
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

# Load 1D diffusion flame
fn_CD = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/CounterDiffusion1D/a=1.70E+04.csv"
f_CD = pd.read_csv(fn_CD, index_col=False)
df_CD = pd.read_csv(fn_CD)
gas_CD = ct.Solution(mech)
fstate_CD = ct.SolutionArray(gas_CD)
fstate_CD.from_pandas(df_CD)
gas_CD.X = "O2:0.21, N2:0.79"; Yox = gas_CD.Y
gas_CD.X = "H2:1.0"; Yfu = gas_CD.Y
zer = mf(gas_CD, Yox, Yfu)
mf_CD = []
for ix, x in enumerate(fstate_CD.grid):
  gas_CD.TPY = fstate_CD.T[ix], fstate_CD.P[ix], fstate_CD.Y[ix,:]
  mf_CD.append(zer.spec2mf(gas_CD.Y))
mf_CD = np.array(mf_CD)

NO_ID = gas_CD.species_index("NO")
NNH_ID = gas_CD.species_index("NNH")
N2O_ID = gas_CD.species_index("N2O")
H2_ID = gas_CD.species_index("H2")
O2_ID = gas_CD.species_index("O2")

# Progress variable definition
coeff_pv = {}
idef = 1
if idef == 0:
  for isp, spn in enumerate(gas1D.species_names):
    if ("N" in spn or spn=="H2" or spn=="O2" or spn=="N2"):
      coeff_pv[spn] = 0.0
    else:
      coeff_pv[spn] = 1.0
elif idef == 1:
  for isp, spn in enumerate(gas1D.species_names):
    if (spn == "H2O"):
      coeff_pv[spn] = 1.0
    elif (spn == "H2"):
      coeff_pv[spn] = -1.0
    elif (spn == "O2"):
      coeff_pv[spn] = -1.0
    else:
      coeff_pv[spn] = 0.0

def get_pveq(zs, equilibrate):
  gas_f = ct.Solution(mech)
  gas_o = ct.Solution(mech)
  gas_m = ct.Solution(mech)
  states = ct.SolutionArray(gas_m)
  print("Input zs for get_pveq: ", zs)
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

  pveq = np.zeros_like(states.T)
  for isp, spn in enumerate(gas_f.species_names):
    pveq = pveq + coeff_pv[spn] * states.Y[:,isp]
  return pveq
pveqs0 = get_pveq(Zs1D, equilibrate = False)
pveqs = get_pveq(Zs1D, equilibrate = True)

#%%
# xyzZC
wt_sum      = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
rho_wtsum   = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
rhoT_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
rhoT2_wtsum = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
mf_wtsum    = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
pv_wtsum    = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
hrr_wtsum   = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
rhoY_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
rhoY2_wtsum = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
wdot_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
wdot2_wtsum = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
#%%
for ifn, fn in enumerate(fns):
  f = h5py.File(fn, 'r+')
  print(fn)
  wt_sum  = wt_sum        + f["DATA"]["volume_mean"]
  rho_wtsum = rho_wtsum   + f["DATA"]["rho_mean"]
  rhoT_wtsum = rhoT_wtsum + f["DATA"]["rhoT_mean"]
  rhoT2_wtsum = rhoT2_wtsum + f["DATA"]["rhoT2_mean"]
  mf_wtsum = mf_wtsum     + f["DATA"]["mixture_fraction_mean"]
  pv_wtsum = pv_wtsum     + f["DATA"]["pv_mean"]
  hrr_wtsum = hrr_wtsum   + f["DATA"]["HeatRelease_mean"]
  for isp in range(0, Nsp):
    fn = "rhoY(" + str(species_names[isp]) + ")_mean"
    rhoY_wtsum[:,:,:,:,:,:,:,isp] = rhoY_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "rhoY2(" + str(species_names[isp]) + ")_mean"
    rhoY2_wtsum[:,:,:,:,:,:,:,isp] = rhoY2_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "wdot(" + str(species_names[isp]) + ")_mean"
    wdot_wtsum[:,:,:,:,:,:,:,isp] = wdot_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "wdot2(" + str(species_names[isp]) + ")_mean"
    wdot2_wtsum[:,:,:,:,:,:,:,isp] = wdot2_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
#%%
# 0 - x
# 1 - y
# 2 - z
# 3 - Z
# 4 - C
# 5 - zone
# 6 - FI
# <|x,zone>
axis_sum = (1, 2, 3, 4)
wt_xzoneFI_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_xzoneFI_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_xzoneFI_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
pv_xzoneFI_wtsum     = np.sum(pv_wtsum, axis=axis_sum)
rhoT_xzoneFI_wtsum   = np.sum(rhoT_wtsum, axis=axis_sum)
rhoT2_xzoneFI_wtsum  = np.sum(rhoT2_wtsum, axis=axis_sum)
hrr_xzoneFI_wtsum    = np.sum(hrr_wtsum, axis=axis_sum)
rhoY_xzoneFI_wtsum   = np.sum(rhoY_wtsum, axis=axis_sum)
rhoY2_xzoneFI_wtsum   = np.sum(rhoY2_wtsum, axis=axis_sum)
wdot_xzoneFI_wtsum   = np.sum(wdot_wtsum, axis=axis_sum)
wdot2_xzoneFI_wtsum   = np.sum(wdot2_wtsum, axis=axis_sum)

#%%
rho_xzoneFI    = rho_xzoneFI_wtsum / wt_xzoneFI_wtsum
mf_xzoneFI     = mf_xzoneFI_wtsum / wt_xzoneFI_wtsum
pv_xzoneFI     = pv_xzoneFI_wtsum / wt_xzoneFI_wtsum
rhoT_xzoneFI   = rhoT_xzoneFI_wtsum / wt_xzoneFI_wtsum
rhoT2_xzoneFI  = rhoT2_xzoneFI_wtsum / wt_xzoneFI_wtsum
hrr_xzoneFI    = hrr_xzoneFI_wtsum / wt_xzoneFI_wtsum
rhoY_xzoneFI   = np.zeros_like(rhoY_xzoneFI_wtsum)
rhoY2_xzoneFI  = np.zeros_like(rhoY2_xzoneFI_wtsum)
wdot_xzoneFI   = np.zeros_like(wdot_xzoneFI_wtsum)
wdot2_xzoneFI  = np.zeros_like(wdot2_xzoneFI_wtsum)
for isp in range(0, Nsp):
  rhoY_xzoneFI[:,:,:,isp] = rhoY_xzoneFI_wtsum[:,:,:,isp] / wt_xzoneFI_wtsum
  rhoY2_xzoneFI[:,:,:,isp] = rhoY2_xzoneFI_wtsum[:,:,:,isp] / wt_xzoneFI_wtsum
  wdot_xzoneFI[:,:,:,isp] = wdot_xzoneFI_wtsum[:,:,:,isp] / wt_xzoneFI_wtsum
  wdot2_xzoneFI[:,:,:,isp] = wdot2_xzoneFI_wtsum[:,:,:,isp] / wt_xzoneFI_wtsum

#%%
x1D = np.linspace(xmin, xmax, nx) + (xmin-xmax)/nx/2
x1D = x1D / 4.5E-4
fig, ax = plt.subplots(figsize=(5,4))
phi = hrr_xzoneFI_wtsum[:,:,:]
phi_sum = np.sum(phi) * (1.28E-2/phi.shape[0])
stitle = r"$\int \mathrm{HRR} \mathrm{d}V_{\mathrm{bin}(x)} / \int \mathrm{HRR} \mathrm{d}V_\mathrm{total}$"
stitle = r"$f_{\mathrm{HRR}}(x \; ; \; \mathcal{M,Z})$"
ax.plot(x1D[:], phi[:,0,0]/phi_sum, color="blue",  linestyle="-", label = r"$\mathrm{NP-Leeward}$")
ax.plot(x1D[:], phi[:,1,0]/phi_sum, color="red",   linestyle="-", label = r"$\mathrm{NP-Windward}$")
ax.plot(x1D[:], phi[:,2,0]/phi_sum, color="green", linestyle="-", label = r"$\mathrm{NP-Interaction}$")
ax.plot(x1D[:], phi[:,0,1]/phi_sum, color="blue",  linestyle="--", label = r"$\mathrm{P-Leeward}$")
ax.plot(x1D[:], phi[:,1,1]/phi_sum, color="red",   linestyle="--", label = r"$\mathrm{P-Windward}$")
ax.plot(x1D[:], phi[:,2,1]/phi_sum, color="green", linestyle="--", label = r"$\mathrm{P-Interaction}$")
ax.set_title(stitle, fontsize=20)
ax.set_xlim([-1, x1D[-1]])
ax.set_xlabel(r"$x/D_j$", fontsize=18)
#ax.set_ylabel(r"$f_{\mathrm{HRR}}(x \; ; \; M,Z)$", fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)
#ax.legend(fontsize=16, bbox_to_anchor=(0.9, 0.45))
plt.savefig("HRR_x.png", dpi=300, bbox_inches="tight")
#%%
lw = 2.0
fig, ax = plt.subplots(figsize=(5,4))
phi = rhoY_xzoneFI_wtsum[:,:,:,NO_ID]
phi_sum = np.sum(phi) * (1.28E-2/phi.shape[0])
stitle = r"$\int \rho Y_\mathrm{NO} \mathrm{d}V_{\mathrm{bin}(x)} / \int \rho Y_\mathrm{NO} \mathrm{d}V_\mathrm{total}$"
stitle = r"$f_{\rho Y_\mathrm{NO}}(x \; ; \; \mathcal{M,R})$"
ax.plot(x1D[:], phi[:,0,0]/phi_sum, color="blue",  linestyle="-", label = r"$\mathrm{NP-Leeward}$", linewidth=lw)
ax.plot(x1D[:], phi[:,1,0]/phi_sum, color="red",   linestyle="-", label = r"$\mathrm{NP-Windward}$", linewidth=lw)
ax.plot(x1D[:], phi[:,2,0]/phi_sum, color="green", linestyle="-", label = r"$\mathrm{NP-Interaction}$", linewidth=lw)
ax.plot(x1D[:], phi[:,0,1]/phi_sum, color="blue",  linestyle="--", label = r"$\mathrm{P-Leeward}$", linewidth=lw)
ax.plot(x1D[:], phi[:,1,1]/phi_sum, color="red",   linestyle="--", label = r"$\mathrm{P-Windward}$", linewidth=lw)
ax.plot(x1D[:], phi[:,2,1]/phi_sum, color="green", linestyle="--", label = r"$\mathrm{P-Interaction}$", linewidth=lw)
ax.set_xlim([-1, x1D[-1]])
ax.set_xlabel(r"$x/D_\mathrm{j}$", fontsize=18)
ax.set_title(stitle,fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)
ax.legend(fontsize=14, frameon = False)
fig.savefig("rhoYNO_x.png", dpi=300)

#%%
fig, ax = plt.subplots(figsize=(5,4))
phi = wdot_xzoneFI_wtsum[:,:,:,NO_ID]
phi_sum = np.sum(wdot_xzoneFI_wtsum[:,:,:,NO_ID]) * (1.28E-2/phi.shape[0])
stitle = r"$\int \rho \dot{\omega}_\mathrm{NO} \mathrm{d}V_{\mathrm{bin}(x)} / \int \rho \dot{\omega}_\mathrm{NO} \mathrm{d}V_{\mathrm{total}}$"
stitle = r"$f_{\rho \dot{\omega}_\mathrm{NO}}(x \; ; \; \mathcal{M,R})$"
ax.plot(x1D[:], phi[:,0,0]/phi_sum, color="blue",  linestyle="-", label = r"$\mathrm{NP-Leeward}$")
ax.plot(x1D[:], phi[:,1,0]/phi_sum, color="red",   linestyle="-", label = r"$\mathrm{NP-Windward}$")
ax.plot(x1D[:], phi[:,2,0]/phi_sum, color="green", linestyle="-", label = r"$\mathrm{NP-Interaction}$")
ax.plot(x1D[:], phi[:,0,1]/phi_sum, color="blue",  linestyle="--", label = r"$\mathrm{P-Leeward}$")
ax.plot(x1D[:], phi[:,1,1]/phi_sum, color="red",   linestyle="--", label = r"$\mathrm{P-Windward}$")
ax.plot(x1D[:], phi[:,2,1]/phi_sum, color="green", linestyle="--", label = r"$\mathrm{P-Interaction}$")
ax.set_xlim([-1, x1D[-1]])
ax.set_xlabel(r"$x/D_\mathrm{j}$", fontsize=18)
ax.set_title(stitle, fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)

#ax.legend()
fig.savefig("rhowdotNO.png", dpi=300)

#%%
fig, ax = plt.subplots(figsize=(5,4))
phi = wdot_xzoneFI_wtsum[:,:,:,NO_ID] / rho_xzoneFI_wtsum[:,:,:] / hrr_xzoneFI_wtsum[:,:,:]
phi_sum = 1.0 #np.sum(phi)
stitle = r"$\int \dot{\omega}_\mathrm{NO} \mathrm{d}V_{\mathrm{bin}(x)} / \int \mathrm{HRR} \mathrm{d}V_{\mathrm{bin}(x)} \; \mathrm{[kg/J]}$"
ax.plot(x1D[:], phi[:,0,0]/phi_sum, color="blue",  linestyle="-", label = r"$\mathrm{NP-Leeward}$")
ax.plot(x1D[:], phi[:,1,0]/phi_sum, color="red",   linestyle="-", label = r"$\mathrm{NP-Windward}$")
ax.plot(x1D[:], phi[:,2,0]/phi_sum, color="green", linestyle="-", label = r"$\mathrm{NP-Interaction}$")
ax.plot(x1D[:], phi[:,0,1]/phi_sum, color="blue",  linestyle="--", label = r"$\mathrm{P-Leeward}$")
ax.plot(x1D[:], phi[:,1,1]/phi_sum, color="red",   linestyle="--", label = r"$\mathrm{P-Windward}$")
ax.plot(x1D[:], phi[:,2,1]/phi_sum, color="green", linestyle="--", label = r"$\mathrm{P-Interaction}$")
ax.set_xlim([0, x1D[-1]])
ax.set_title(stitle)
ax.legend()
fig.savefig("wdot2hrr_x.png", dpi=300)

#%%
fig, ax = plt.subplots(figsize=(5,4))
phi = -wdot_xzoneFI_wtsum[:,:,:,NO_ID] / wdot_xzoneFI_wtsum[:,:,:,H2_ID]
phi_sum = 1.0 #np.sum(phi)
stitle = r"$-\int \dot{\omega}_\mathrm{NO} \mathrm{d}V_{\mathrm{bin}(x)} / \int \dot{\omega}_\mathrm{H2} \mathrm{d}V_{\mathrm{bin}(x)} \; \mathrm{[kg/kg]}$"
ax.plot(x1D[:], phi[:,0,0]/phi_sum, color="blue",  linestyle="-", label = r"$\mathrm{NP-Leeward}$")
ax.plot(x1D[:], phi[:,1,0]/phi_sum, color="red",   linestyle="-", label = r"$\mathrm{NP-Windward}$")
ax.plot(x1D[:], phi[:,2,0]/phi_sum, color="green", linestyle="-", label = r"$\mathrm{NP-Interaction}$")
ax.plot(x1D[:], phi[:,0,1]/phi_sum, color="blue",  linestyle="--", label = r"$\mathrm{P-Leeward}$")
ax.plot(x1D[:], phi[:,1,1]/phi_sum, color="red",   linestyle="--", label = r"$\mathrm{P-Windward}$")
ax.plot(x1D[:], phi[:,2,1]/phi_sum, color="green", linestyle="--", label = r"$\mathrm{P-Interaction}$")
ax.set_xlim([0, x1D[-1]])
ax.set_title(stitle)
ax.legend()
fig.savefig("testlatex.png", dpi=300)

#%%
