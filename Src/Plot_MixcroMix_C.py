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
field_names = [["T",            "Y(H2)",    "Y(H)",    "Y(H2O2)",     "Y(NO)"],
                ["HeatRelease", "wdot(H2)", "wdot(H)", "wdot(H2O2)",  "wdot(NO)"],
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
# <|Z>
axis_sum = (0, 1, 2, 3, 5, 6)
wt_Z_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_Z_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_Z_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
# <|Z,C>
axis_sum = (0, 1, 2, 5, 6)
wt_ZC_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_ZC_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_ZC_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
pv_ZC_wtsum     = np.sum(pv_wtsum, axis=axis_sum)
rhoT_ZC_wtsum   = np.sum(rhoT_wtsum, axis=axis_sum)
rhoT2_ZC_wtsum  = np.sum(rhoT2_wtsum, axis=axis_sum)
hrr_ZC_wtsum    = np.sum(hrr_wtsum, axis=axis_sum)
hrr2_ZC_wtsum   = np.sum(hrr2_wtsum, axis=axis_sum)
rhoY_ZC_wtsum   = np.sum(rhoY_wtsum, axis=axis_sum)
rhoY2_ZC_wtsum   = np.sum(rhoY2_wtsum, axis=axis_sum)
wdot_ZC_wtsum   = np.sum(wdot_wtsum, axis=axis_sum)
wdot2_ZC_wtsum   = np.sum(wdot2_wtsum, axis=axis_sum)

#%%
rho_ZC    = rho_ZC_wtsum / wt_ZC_wtsum
mf_ZC     = mf_ZC_wtsum / wt_ZC_wtsum
pv_ZC     = pv_ZC_wtsum / wt_ZC_wtsum
rhoT_ZC   = rhoT_ZC_wtsum / wt_ZC_wtsum
rhoT2_ZC  = rhoT2_ZC_wtsum / wt_ZC_wtsum
hrr_ZC    = hrr_ZC_wtsum / wt_ZC_wtsum
hrr2_ZC   = hrr2_ZC_wtsum / wt_ZC_wtsum
rhoY_ZC   = np.zeros_like(rhoY_ZC_wtsum)
rhoY2_ZC  = np.zeros_like(rhoY2_ZC_wtsum)
wdot_ZC   = np.zeros_like(wdot_ZC_wtsum)
wdot2_ZC  = np.zeros_like(wdot2_ZC_wtsum)
for isp in range(0, Nsp):
  rhoY_ZC[:,:,isp] = rhoY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  rhoY2_ZC[:,:,isp] = rhoY2_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  wdot_ZC[:,:,isp] = wdot_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  wdot2_ZC[:,:,isp] = wdot2_ZC_wtsum[:,:,isp] / wt_ZC_wtsum

fig, ax = plt.subplots(figsize=(3.6,4))
phi = pv_ZC
ax.imshow(phi.transpose(), origin="lower",
          vmin=0, vmax=1,
          cmap = "jet",
          extent = [0.0, 1.0, 0.0, 1.2])
ax.set_xlabel(r"$Z$", fontsize=labelsize)
ax.set_xlim([0.0, 1.0])
ax.set_ylabel(r"$C$", fontsize=labelsize)
ax.set_title(r"$\langle T \; | \; Z, C \rangle~\mathrm{[K]}$", fontsize=labelsize+6)
#%%
Zt = 0.025
iz1D = np.argmax((Zs1D > Zt))
iz3D = np.argmax((zout > Zt))-1
fn1D = fns1D_sorted[iz1D]
df = pd.read_csv(fn1D)
fstate = ct.SolutionArray(gas1D)
fstate.from_pandas(df)
net_production_rates1D = fstate.net_production_rates
partial_molar_enthalpies1D = fstate.partial_molar_enthalpies
hrr1D = -np.sum(net_production_rates1D * partial_molar_enthalpies1D, axis=(1))

pv1D = np.zeros_like(fstate.T)
for isp, spn in enumerate(gas1D.species_names):
  pv1D = pv1D + coeff_pv[spn] * fstate.Y[:, isp]
pv1D = (pv1D - pveqs0[iz1D]) / (pveqs[iz1D] - pveqs0[iz1D])

# Plot parameter
npy = len(field_names); npx = len(field_names[0])
fig_unit_y = 2
labelsize = 16
figsize = (fig_unit_y*npx*1.45, fig_unit_y*npy*1.05)
fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)
for ipy in range(0, npy):
  for ipx in range(0, npx):
    ax = axs[ipy, ipx]
    field_name = field_names[ipy][ipx]
    pv3D = pv_ZC[iz3D,:]
    if (field_name == "T"):
      # 1D
      ax.plot(pv1D, fstate.T, color="r", linestyle="--", linewidth=lw, label = r"$\mathrm{1D}$")
      # 3D
      T_favm = rhoT_ZC[iz3D,:] / rho_ZC[iz3D,:]
      ax.plot(pv3D, T_favm, color="k", linestyle="-", linewidth=lw, label=r"$\langle T | C \rangle$")
      # 3D - std
      rhoT_fms = rhoT2_ZC[iz3D,:] - rho_ZC[iz3D,:]*T_favm*T_favm
      T_fms = np.sqrt(rhoT_fms / rho_ZC[iz3D,:])
      y1 = T_favm - T_fms
      y2 = T_favm + T_fms
      ax.fill_between(pv3D, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$T\;\mathrm{[K]}$", fontsize=labelsize)
      ax.set_ylim([500, 2500])
      #ax.legend(fontsize=labelsize)
    elif (field_name[0:2] == "Y("):
      spn = re.split("\(", field_name)[-1]
      spn = re.split("\)", spn)[0]
      ispp = gas1D.species_index(spn)
      if spn == "H2":
        spl = "H_2"
      elif spn == "H2O2":
        spl = "H_2 O_2"
      else:
        spl = spn
      # 1D
      Y_1D   = fstate.Y[:,ispp]
      ax.plot(pv1D, Y_1D, color="r", linestyle="--", linewidth=lw)
      # 3D
      Y_favm = rhoY_ZC[iz3D,:,ispp] / rho_ZC[iz3D,:]
      ax.plot(pv3D, Y_favm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      rhoY_fms = rhoY2_ZC[iz3D,:,ispp] - rho_ZC[iz3D,:]*Y_favm*Y_favm
      Y_fms = np.sqrt(rhoY_fms / rho_ZC[iz3D,:])
      y1 = Y_favm - Y_fms
      y2 = Y_favm + Y_fms
      ax.fill_between(pv3D, y1, y2, color="gray", alpha=0.5)
      indx = (Y_favm<1) & (Y_favm>0)
      Y_max_1D = np.amax(Y_1D)
      Y_max_3D = np.amax(Y_favm[indx])
      Yrms_max_3D = np.amax(Y_fms[indx])
      ax.set_ylim([0, np.amax([Y_max_3D+Yrms_max_3D/2, Y_max_3D])*1.15])
      ax.set_title(r"$Y_\mathrm{"+spl+"}$", fontsize=labelsize)
      if spn == "H2":
        ax.set_ylim([0, 0.12])
    elif (field_name[0:4] == "wdot"):
      spn = re.split("\(", field_name)[-1]
      spn = re.split("\)", spn)[0]
      ispp = gas1D.species_index(spn)
      if spn == "H2":
        spl = "H_2"
      elif spn == "H2O2":
        spl = "H_2 O_2"
      else:
        spl = spn
      # 1D
      wdot1D = net_production_rates1D[:,ispp] * gas1D.molecular_weights[ispp]
      ax.plot(pv1D, wdot1D, color="r", linestyle="--", linewidth=lw)
      # 3D
      wdot_favm = wdot_ZC[iz3D,:,ispp] / rho_ZC[iz3D,:]
      ax.plot(pv3D, wdot_favm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      wdot_fms = wdot2_ZC[iz3D,:,ispp] - rho_ZC[iz3D,:]*wdot_favm*wdot_favm
      wdot_fms = np.sqrt(wdot_fms)
      y1 = wdot_favm + wdot_fms
      y2 = wdot_favm - wdot_fms
      ax.fill_between(pv3D, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$\dot{\omega}_\mathrm{"+spl+"}\; \mathrm{[kg/m^3 s]}$", fontsize=labelsize)
    elif (field_name[0:4] == "Heat"):
      # 1D
      ax.plot(pv1D, hrr1D, color="r", linestyle="--", linewidth=lw)
      # 3D
      #wdot_favm = wdot_ZC[iz3D,:,ispp] / rho_ZC[iz3D,:]
      #ax.plot(pv3D, wdot_favm, color="k", linestyle="-")
      hrr_ensm = hrr_ZC[iz3D,:]
      ax.plot(pv3D, hrr_ensm, color="k", linestyle="-", linewidth=lw)
      # 3D - std
      hrr2_ensm = hrr2_ZC[iz3D,:]
      hrr_rms = np.sqrt(hrr2_ensm - hrr_ensm*hrr_ensm)
      y1 = hrr_ensm - hrr_rms
      y2 = hrr_ensm + hrr_rms
      ax.fill_between(pv3D, y1, y2, color="gray", alpha=0.5)
      ax.set_title(r"$\mathrm{HRR}\; \mathrm{[J/m^3 s]}$", fontsize=labelsize)

    if (ipy == npy-1):
      ax.set_xlabel(r"$C$", fontsize = labelsize)
      ax.set_xticks(np.array([0, 0.5, 1]))
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

plt.savefig("FlameStructure_C_Z=" + "%.3E"%(Zt)+".png", dpi=300, bbox_inches="tight")
#%%