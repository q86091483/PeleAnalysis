#%%
import h5py
import numpy as np
from numpy import genfromtxt
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
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix_RPA/"
# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_1*.h5")
#%%
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

Nr = 235

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
spn = "NO"
NO_ID = gas1D.species_index(spn)
spn = "N2"
N2_ID = gas1D.species_index(spn)
spn = "N"
N_ID = gas1D.species_index(spn)
spn = "N2O"
N2O_ID = gas1D.species_index(spn)
spn = "NO2"
NO2_ID = gas1D.species_index("NO2")
spn = "NNH"
NNH_ID = gas1D.species_index("NNH")
spn = "NH"
NH_ID = gas1D.species_index("NH")
spn = "NH2"
NH2_ID = gas1D.species_index("NH2")
spn = "HNO"
HNO_ID = gas1D.species_index("HNO")
spn = "HON"
HON_ID = gas1D.species_index("HON")
spn = "HONO"
HONO_ID = gas1D.species_index("HONO")
spn = "HNO2"
HNO2_ID = gas1D.species_index("HNO2")
spn = "HNO3"
HNO3_ID = gas1D.species_index("HNO3")

Nsp = len(gas1D.species_names)
species_names = gas1D.species_names

# Load transfer csv
csv_path = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/nuig_H2_4atm/n-N-sp1-sp2_nuig_H2_32sp_4atm.csv"
transfer_matrix = genfromtxt(csv_path, delimiter=',')
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
prr_wtsum   = np.zeros((nx, ny, nz, nmf, npv, nzo, nfi, Nr))

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
  for ir in range(0, Nr):
    fn = "R(" + str(ir) + ")_mean"
    prr_wtsum[:,:,:,:,:,:,:,ir] = prr_wtsum[:,:,:,:,:,:,:,ir] + f["DATA"][fn] 
#%%
# 0 - x
# 1 - y
# 2 - z
# 3 - Z
# 4 - C
# 5 - zone
# 6 - FI
# <|Z,C>
axis_sum = (1, 2, 3, 4)
wt_xzoneFI_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_xzoneFI_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_xzoneFI_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
pv_xzoneFI_wtsum     = np.sum(pv_wtsum, axis=axis_sum)
rhoT_xzoneFI_wtsum   = np.sum(rhoT_wtsum, axis=axis_sum)
hrr_xzoneFI_wtsum    = np.sum(hrr_wtsum, axis=axis_sum)
rhoY_xzoneFI_wtsum   = np.sum(rhoY_wtsum, axis=axis_sum)
rhoY2_xzoneFI_wtsum   = np.sum(rhoY2_wtsum, axis=axis_sum)
wdot_xzoneFI_wtsum   = np.sum(wdot_wtsum, axis=axis_sum)
wdot2_xzoneFI_wtsum   = np.sum(wdot2_wtsum, axis=axis_sum)
prr_xzoneFI_wtsum   = np.sum(prr_wtsum, axis=axis_sum)

# progressRateVolume[0:NUM_REACTIONS]
progressRateVolume_domainSum = np.sum(prr_wtsum, axis=(0,1,2,3,4,5,6))
# Domain volume
Volume_domain = np.sum(wt_sum, axis=(0,1,2,3,4,5,6))
#%%
prr = progressRateVolume_domainSum / Volume_domain
ix = 4
iy = 0
iz = 0
iZ = 1
iC = 8
izone = 0
iFI = 0
prr = prr_wtsum[ix, iy, iz, iZ, iC, izone, iFI,:] / wt_sum[ix, iy, iz, iZ, iC, izone, iFI]

#%% Whole domain
prr = np.sum(prr_wtsum, axis=(0,1,2,3,4,5,6)) # prr[nx,ny,nz,nZ,nC,nzone,nFI,NUM_REACTIONs]
volume = np.sum(wt_sum, axis=(0,1,2,3,4,5,6)) #
graph = np.zeros((Nsp, Nsp))
for iline in range(0, transfer_matrix.shape[0]):
  ir    = int(transfer_matrix[iline, 0])
  isp0  = int(transfer_matrix[iline, 1]) 
  isp1  = int(transfer_matrix[iline, 2])
  natom = int(transfer_matrix[iline, 3])
  graph[isp0, isp1] = graph[isp0, isp1] + prr[ir] * natom
graph = graph - graph.T
#graph[graph<0] = 0
target_ID = NH_ID
target_name = "NH"
for isp in range(0, Nsp):
  print(target_name + " to ".rjust(5)+gas1D.species_names[isp].ljust(5), ": " 
         "%.2E" % (graph[target_ID, isp] / np.sum(graph[N2_ID, :])))
for isp in range(0, Nsp):
  print(target_name + " from ".rjust(5)+gas1D.species_names[isp].ljust(5), ": " 
         "%.2E" % (graph[isp, target_ID] / np.sum(graph[N2_ID, :])))



#%% Premixed
prr_mode = np.sum(prr_wtsum, axis=(0,1,2,3,4,5)) # prr[nx,ny,nz,nZ,nC,nzone,nFI,NUM_REACTIONs]
volume_mode = np.sum(wt_sum, axis=(0,1,2,3,4,5)) #
imode = 1 # Premixed
graph = np.zeros((Nsp, Nsp))
for iline in range(0, transfer_matrix.shape[0]):
  ir    = int(transfer_matrix[iline, 0])
  isp0  = int(transfer_matrix[iline, 1]) 
  isp1  = int(transfer_matrix[iline, 2])
  natom = int(transfer_matrix[iline, 3])
  graph[isp0, isp1] = graph[isp0, isp1] + prr_mode[imode,ir] * natom
graph = graph - graph.T
graph[graph<0] = 0
target_name = "N2"
target_ID = gas1D.species_index(target_name)
for isp in range(0, Nsp):
  print(target_name + " to ".rjust(5)+gas1D.species_names[isp].ljust(5), ": " 
         "%.2E" % (graph[target_ID, isp] / np.sum(graph[N2_ID, :])))
for isp in range(0, Nsp):
  print(target_name + " from ".rjust(5)+gas1D.species_names[isp].ljust(5), ": " 
         "%.2E" % (graph[isp, target_ID] / np.sum(graph[N2_ID, :])))

#%%
fig, ax = plt.subplots(figsize=(5, 4))
spn_rpa = ["NO", "N", "NNH", "N_2O"]
ixs = [1, 2, 3, 4]
ratio_premixed = [5.72E-2, 5.78E-2, 1E-2, 8.84E-1]
ratio_nonpremixed = [3.99E-1, 3.97E-1, 1.56E-1, 4.51E-2]

fig, ax = plt.subplots()
ax.plot(ixs, ratio_nonpremixed, color="r", linewidth=3.0, 
        marker="o", markersize=15.0, label=r"$\mathrm{Nonpremixed}$" )
ax.plot(ixs, ratio_premixed, color="b", linewidth=3.0, 
        marker="o", markersize=15.0, label=r"$\mathrm{Premixed}$")
ax.legend(fontsize=18)
xticklabels = []
for i in range(0, len(spn_rpa)):
  xticklabels.append(r"$\mathrm{"+spn_rpa[i]+"}$")
ax.set_xticks(ixs)
ax.set_xticklabels(xticklabels)
ax.set_xlabel(r"$\mathrm{Species}$", fontsize = 18)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=16)
ax.set_ylim([0,1])
ax.set_ylabel("Ratio of N flux",
              fontsize = 20)

#%%
bar_width = 0.1
multiplier = 0
for ibar, spn in enumerate(spn_rpa):
  offset = bar_width*ibar + bar_width * multiplier
  rects = ax.bar(offset, ratio_nonpremixed[ibar], bar_width, label=spn_rpa[ibar])
  #rects = ax.bar(offset*2, ratio_premixed[ibar], bar_width, label=spn_rpa[ibar])

#%%

#%%
import graphviz
f = graphviz.Digraph()

spo = "N2"; i = gas1D.species_index(spo)
origin = f"{spo}"
spd = "NO"; j = gas1D.species_index(spd)
destin = f"{spd}"
color = "black"
width = graph[i,j] / np.sum(graph[i,:])
#width = max(1, width)
#label = f"  {100*graph[i, j]/graph.max():.2f} %\n\n"
label = f"   {width:.2e}\n\n"
f.edge(origin, destin, label=label, penwidth=str(width), color=color)
f

# %%
