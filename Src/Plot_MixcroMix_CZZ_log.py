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
from mixture_fraction import mf

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Input
# Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix_ZC/"
# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_1*.h5")
# Data folder - where 1D flame data is stored
data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/UnstrainedPremixed1D"
# Field names
field_names = [["T",            "Y(H2)",    "Y(H)",    "Y(NO)",     "Y(NNH)", "Y(N2O)"], 
                ["HeatRelease", "wdot(H2)", "wdot(H)", "wdot(NO)",  "wdot(NNH)", "wdot(N2O)"],
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
labelsize = 24
lw = 3
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
hrr2_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi))
rhoY_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
rhoY2_wtsum = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
wdot_wtsum  = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))
wdot2_wtsum = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nsp))

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
  hrr2_wtsum = hrr2_wtsum   + f["DATA"]["HeatRelease2_mean"]
  for isp in range(0, Nsp):
    fn = "rhoY(" + str(species_names[isp]) + ")_mean"
    rhoY_wtsum[:,:,:,:,:,:,:,isp] = rhoY_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "rhoY2(" + str(species_names[isp]) + ")_mean"
    rhoY2_wtsum[:,:,:,:,:,:,:,isp] = rhoY2_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "wdot(" + str(species_names[isp]) + ")_mean"
    wdot_wtsum[:,:,:,:,:,:,:,isp] = wdot_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]
    fn = "wdot2(" + str(species_names[isp]) + ")_mean"
    wdot2_wtsum[:,:,:,:,:,:,:,isp] = wdot2_wtsum[:,:,:,:,:,:,:,isp] + f["DATA"][fn]

#%% <|Z>
axis_sum = (0, 1, 2, 4, 5, 6)
wt_Z_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_Z_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_Z_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
pv_Z_wtsum     = np.sum(pv_wtsum, axis=axis_sum)
rhoT_Z_wtsum   = np.sum(rhoT_wtsum, axis=axis_sum)
rhoT2_Z_wtsum  = np.sum(rhoT2_wtsum, axis=axis_sum)
hrr_Z_wtsum    = np.sum(hrr_wtsum, axis=axis_sum)
hrr2_Z_wtsum   = np.sum(hrr2_wtsum, axis=axis_sum)
rhoY_Z_wtsum   = np.sum(rhoY_wtsum, axis=axis_sum)
rhoY2_Z_wtsum   = np.sum(rhoY2_wtsum, axis=axis_sum)
wdot_Z_wtsum   = np.sum(wdot_wtsum, axis=axis_sum)
wdot2_Z_wtsum   = np.sum(wdot2_wtsum, axis=axis_sum)

rho_Z    = rho_Z_wtsum / wt_Z_wtsum 
mf_Z     = mf_Z_wtsum / wt_Z_wtsum 
pv_Z     = pv_Z_wtsum / wt_Z_wtsum 
rhoT_Z   = rhoT_Z_wtsum / wt_Z_wtsum 
rhoT2_Z  = rhoT2_Z_wtsum / wt_Z_wtsum 
hrr_Z    = hrr_Z_wtsum / wt_Z_wtsum 
hrr2_Z   = hrr2_Z_wtsum / wt_Z_wtsum 
rhoY_Z   = np.zeros_like(rhoY_Z_wtsum)
rhoY2_Z  = np.zeros_like(rhoY2_Z_wtsum)
wdot_Z   = np.zeros_like(wdot_Z_wtsum)
wdot2_Z  = np.zeros_like(wdot2_Z_wtsum)
for isp in range(0, Nsp):
  rhoY_Z[:,isp] = rhoY_Z_wtsum[:,isp] / wt_Z_wtsum
  rhoY2_Z[:,isp] = rhoY2_Z_wtsum[:,isp] / wt_Z_wtsum
  wdot_Z[:,isp] = wdot_Z_wtsum[:,isp] / wt_Z_wtsum
  wdot2_Z[:,isp] = wdot2_Z_wtsum[:,isp] / wt_Z_wtsum

# Load 1D diffusion flame - 1
#%%
fn_CD = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/CounterDiffusion1D/a=1.44E+04.csv"
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
#%%

net_production_rates1D = fstate_CD.net_production_rates
partial_molar_enthalpies1D = fstate_CD.partial_molar_enthalpies
hrr1D = -np.sum(net_production_rates1D * partial_molar_enthalpies1D, axis=(1))

# Load 1D diffusion flame - 2
fn_CD_1 = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/CounterDiffusion1D/a=1.44E+04.csv"
f_CD_1 = pd.read_csv(fn_CD_1, index_col=False)
df_CD_1 = pd.read_csv(fn_CD_1)
gas_CD_1 = ct.Solution(mech)
fstate_CD_1 = ct.SolutionArray(gas_CD_1)
fstate_CD_1.from_pandas(df_CD_1)
gas_CD_1.X = "O2:0.21, N2:0.79"; Yox = gas_CD_1.Y
gas_CD_1.X = "H2:1.0"; Yfu = gas_CD_1.Y
zer = mf(gas_CD_1, Yox, Yfu)
mf_CD_1 = []
for ix, x in enumerate(fstate_CD_1.grid):
  gas_CD_1.TPY = fstate_CD_1.T[ix], fstate_CD_1.P[ix], fstate_CD_1.Y[ix,:]
  mf_CD_1.append(zer.spec2mf(gas_CD_1.Y))
mf_CD_1 = np.array(mf_CD_1) 
net_production_rates1D_1 = fstate_CD_1.net_production_rates
partial_molar_enthalpies1D_1 = fstate_CD_1.partial_molar_enthalpies
hrr1D_1 = -np.sum(net_production_rates1D_1 * partial_molar_enthalpies1D_1, axis=(1))
#%%
# Plot parameter
npy = len(field_names); npx = len(field_names[0])
fig_unit_y = 2
labelsize = 16
figsize = (fig_unit_y*npx*1.35, fig_unit_y*npy*1.05)
fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)
for ipy in range(0, npy):
  for ipx in range(0, npx):
    ax = axs[ipy, ipx]
    field_name = field_names[ipy][ipx]
    if (field_name == "T"):
      # 1D
      ax.plot(mf_CD, fstate_CD.T, color="r", linestyle="--", linewidth=lw, label = r"$\mathrm{1D}$")
      ax.plot(mf_CD_1, fstate_CD_1.T, color="b", linestyle="--", linewidth=lw, label = r"$\mathrm{1D}$")
      # 3D
      T_favm = rhoT_Z[:] / rho_Z[:]
      ax.plot(zout, T_favm, color="k", linestyle="-", linewidth=lw, label=r"$\langle T | C \rangle$") 
      # 3D - std
      rhoT_fms = rhoT2_Z[:] - rho_Z[:]*T_favm*T_favm
      T_fms = np.sqrt(rhoT_fms / rho_Z[:])
      y1 = T_favm - T_fms
      y2 = T_favm + T_fms
      ax.fill_between(zout, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$T\;\mathrm{[K]}$", fontsize=labelsize)
      ax.set_ylim([300, 2800])
      if True:
        ax.set_xscale('log')
      ax.set_xlim([0.01, 1.0])
      ax.legend(fontsize=labelsize)
    elif (field_name[0:2] == "Y("):
      spn = re.split("\(", field_name)[-1]
      spn = re.split("\)", spn)[0]
      ispp = gas1D.species_index(spn)
      # 1D
      Y_1D   = fstate_CD.Y[:,ispp]
      ax.plot(mf_CD, Y_1D, color="r", linestyle="--", linewidth=lw)
      Y_1D_1 = fstate_CD_1.Y[:,ispp]
      ax.plot(mf_CD_1, Y_1D_1, color="b", linestyle="--", linewidth=lw)
      # 3D
      Y_favm = rhoY_Z[:,ispp] / rho_Z[:]
      ax.plot(zout, Y_favm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      rhoY_fms = rhoY2_Z[:,ispp] - rho_Z[:]*Y_favm*Y_favm
      Y_fms = np.sqrt(rhoY_fms / rho_Z[:])
      y1 = Y_favm - Y_fms
      y2 = Y_favm + Y_fms
      ax.fill_between(zout, y1, y2, color="gray", alpha=0.5)
      indx = (Y_favm<1) & (Y_favm>0) 
      Y_max_1D = np.amax(Y_1D)
      Y_max_3D = np.amax(Y_favm[indx])
      Yrms_max_3D = np.amax(Y_fms[indx])
      ax.set_ylim([0, np.amax([Y_max_3D+Yrms_max_3D/2, Y_max_3D])*1.1])
      ax.set_title(r"$Y_\mathrm{"+spn+"}$", fontsize=labelsize)
      if True:
        ax.set_xscale('log')
      ax.set_xlim([0.01, 1.0])
      if spn == "NO":
        ax.set_xlim([0, 1.0])
    elif (field_name[0:4] == "wdot"):
      spn = re.split("\(", field_name)[-1]
      spn = re.split("\)", spn)[0]
      ispp = gas1D.species_index(spn)
      # 1D
      wdot1D = net_production_rates1D[:,ispp] * gas1D.molecular_weights[ispp]
      ax.plot(mf_CD, wdot1D, color="r", linestyle="--", linewidth=lw)
      wdot1D_1 = net_production_rates1D_1[:,ispp] * gas1D.molecular_weights[ispp]
      ax.plot(mf_CD_1, wdot1D_1, color="b", linestyle="--", linewidth=lw)
      # 3D
      wdot_favm = wdot_Z[:,ispp] / rho_Z[:]
      ax.plot(zout, wdot_favm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      wdot_fms = wdot2_Z[:,ispp] - rho_Z[:]*wdot_favm*wdot_favm
      wdot_fms = np.sqrt(wdot_fms)
      y1 = wdot_favm + wdot_fms
      y2 = wdot_favm - wdot_fms
      ax.fill_between(zout, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$\dot{\omega}_\mathrm{"+spn+"}\; \mathrm{[kg/m^3 s]}$", fontsize=labelsize)
      if True:
        ax.set_xscale('log')
      ax.set_xlim([0.01, 1.0])
      ax.set_xticks([0.01, 1.0])
      #if spn == "NO":
      #  ax.set_xlim([0.01, 1.0])
    elif (field_name[0:4] == "Heat"):
      # 1D
      ax.plot(mf_CD, hrr1D, color="r", linestyle="--", linewidth=lw)
      ax.plot(mf_CD_1, hrr1D_1, color="b", linestyle="--", linewidth=lw)
      # 3D
      #wdot_favm = wdot_ZC[iz3D,:,ispp] / rho_ZC[iz3D,:]
      #ax.plot(pv3D, wdot_favm, color="k", linestyle="-")
      hrr_ensm = hrr_Z[:]
      ax.plot(zout, hrr_ensm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      hrr2_ensm = hrr2_Z[:]
      hrr_rms = np.sqrt(hrr2_ensm - hrr_ensm*hrr_ensm)
      y1 = hrr_ensm - hrr_rms
      y2 = hrr_ensm + hrr_rms
      ax.fill_between(zout, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$\mathrm{HRR}\; \mathrm{[J/m^3 s]}$", fontsize=labelsize)
      ax.set_ylim([-5E10, 1E11])
      if True:
        ax.set_xscale('log')
      ax.set_xlim([0.001, 1.0])

    if (ipy == npy-1):
      ax.set_xlabel(r"$Z$", fontsize = labelsize)
      ax.set_xticks(np.array([0.001, 0.01, 0.1, 1.0]))
      ax.set_xticks
    else:
      #ax.set_xlabel([])
      ax.set_xticks(np.array([]))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(3,1))
    tit = ax.yaxis.get_offset_text()
    tit.set_x(-0.15)
    ax.yaxis.get_offset_text().set_fontsize(labelsize)

    #if (ipx == 0):
      #ax.set_ylabel(r"$z / D_{j}$", fontsize = labelsize)
      #ax.set_yticks(np.array([0, 5, 10]))
    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    ax.tick_params(axis='both', which='minor', labelsize=labelsize)
    #ax.set_xlim([0, 1.0])
      
plt.savefig("FlameStructure_C_Z=" + ".png", dpi=300, bbox_inches="tight")
#%%