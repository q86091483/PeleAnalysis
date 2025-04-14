#%%
import h5py
import numpy as np
from numpy import genfromtxt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import glob as glob
import os
import os.path as path
import re
import pandas as pd
import cantera as ct
from mixture_fraction import mf

plt.rcParams['text.usetex'] = True
os.environ["PATH"] += os.pathsep + "/scinet/niagara/software/2019b/opt/base/texlive/20230314/bin/x86_64-linux/"
labelsize = 22; lw = 1.5

# Input - Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_fluxRPA/"
fns = glob.glob(case_folder + "plt_1*.h5")
def get_Z_key(s):
  ss = re.split("/", s)[-1]
  ss = re.split("_", ss)[1]
  return(float(ss))
fns_sorted = sorted(fns, key=get_Z_key)

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
ntres = tet.shape[7]
Nr = 235
zout = np.linspace(0, 1, nmf)
pvout = np.linspace(0, 1, npv)

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
gas1D = ct.Solution(mech,  transport_model="Multi")
Nsp = len(gas1D.species_names)
species_names = gas1D.species_names
spn = "NO";   NO_ID = gas1D.species_index(spn)
spn = "N2";   N2_ID = gas1D.species_index(spn)
spn = "N";    N_ID = gas1D.species_index(spn)
spn = "N2O";  N2O_ID = gas1D.species_index(spn)
spn = "NO2";  NO2_ID = gas1D.species_index(spn)
spn = "NNH";  NNH_ID = gas1D.species_index(spn)
spn = "NH";   NH_ID = gas1D.species_index(spn)
spn = "NH2";  NH2_ID = gas1D.species_index(spn)
spn = "HNO";  HNO_ID = gas1D.species_index(spn)
spn = "HON";  HON_ID = gas1D.species_index(spn)
spn = "HONO"; HONO_ID = gas1D.species_index(spn)
spn = "HNO2"; HNO2_ID = gas1D.species_index(spn)
spn = "HNO3"; HNO3_ID = gas1D.species_index(spn)
spn = "H2";   H2_ID = gas1D.species_index(spn)
spn = "H";    H_ID = gas1D.species_index(spn)
spn = "O2";   O2_ID = gas1D.species_index(spn)
spn = "OH";   OH_ID = gas1D.species_index(spn)
spn = "H2O";  H2O_ID = gas1D.species_index(spn)
spn = "H2O2"; H2O2_ID = gas1D.species_index(spn)
spn = "HO2";  HO2_ID = gas1D.species_index(spn)

#%% Load transfer csv
csv_path = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/nuig_H2_4atm/n-N-sp1-sp2_nuig_H2_32sp_4atm.csv"
transfer_matrix = genfromtxt(csv_path, delimiter=',')
for i, line in enumerate(transfer_matrix):
  ir = int(line[0])
  isp0 = int(line[1])
  isp1 = int(line[2])
  natom = int(line[3])
  print(str(gas1D.reaction(ir)).ljust(25), ", ", line)
  strc = "  " + species_names[isp0] + " -> " + species_names[isp1]
  print(strc.rjust(13), natom)

#%% Load 1D premixed flames
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

#%% Load 0D equilibirum state
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

#%% Load 1D diffusion flame
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

#%% Read x, y, z, Z, C, zone, FI, tres
dim_cond = (nx, ny, nz, nmf, npv, nzo, nfi, ntres)
dim_cond_sp = dim_cond + (Nsp,)
dim_cond_r = dim_cond + (Nr,)
wt_sum      = np.zeros(dim_cond)
rho_wtsum   = np.zeros(dim_cond)
rhoT_wtsum  = np.zeros(dim_cond)
rhoT2_wtsum = np.zeros(dim_cond)
mf_wtsum    = np.zeros(dim_cond)
pv_wtsum    = np.zeros(dim_cond)
hrr_wtsum   = np.zeros(dim_cond)
hrr2_wtsum  = np.zeros(dim_cond)
agepv_wtsum = np.zeros(dim_cond)
rhoY_wtsum  = np.zeros(dim_cond_sp)
rhoY2_wtsum = np.zeros(dim_cond_sp)
wdot_wtsum  = np.zeros(dim_cond_sp)
wdot2_wtsum = np.zeros(dim_cond_sp)
prr_wtsum   = np.zeros(dim_cond_r)
diffY_wtsum = np.zeros(dim_cond_sp)
diffWY_wtsum = np.zeros(dim_cond_sp)
diffTY_wtsum = np.zeros(dim_cond_sp)
convY_wtsum = np.zeros(dim_cond_sp)

for ifn, fn in enumerate(fns_sorted):
  f = h5py.File(fn, 'r+')
  wt_sum  = wt_sum        + f["DATA"]["volume_mean"]
  rho_wtsum = rho_wtsum   + f["DATA"]["rho_mean"]
  rhoT_wtsum = rhoT_wtsum + f["DATA"]["rhoT_mean"]
  rhoT2_wtsum = rhoT2_wtsum + f["DATA"]["rhoT2_mean"]
  mf_wtsum = mf_wtsum     + f["DATA"]["mixture_fraction_mean"]
  pv_wtsum = pv_wtsum     + f["DATA"]["pv_mean"]
  hrr_wtsum = hrr_wtsum   + f["DATA"]["HeatRelease_mean"]
  hrr2_wtsum = hrr2_wtsum   + f["DATA"]["HeatRelease2_mean"]
  agepv_wtsum = agepv_wtsum + f["DATA"]["agepv_1_mean"]
  for isp in range(0, Nsp):
    fdn = "rhoY(" + str(species_names[isp]) + ")_mean"
    rhoY_wtsum[:,:,:,:,:,:,:,:,isp] = rhoY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "rhoY2(" + str(species_names[isp]) + ")_mean"
    rhoY2_wtsum[:,:,:,:,:,:,:,:,isp] = rhoY2_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "wdot(" + str(species_names[isp]) + ")_mean"
    wdot_wtsum[:,:,:,:,:,:,:,:,isp] = wdot_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "wdot2(" + str(species_names[isp]) + ")_mean"
    wdot2_wtsum[:,:,:,:,:,:,:,:,isp] = wdot2_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "diffY(" + str(species_names[isp]) + ")_mean"
    diffY_wtsum[:,:,:,:,:,:,:,:,isp] = diffY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "diffWY(" + str(species_names[isp]) + ")_mean"
    diffWY_wtsum[:,:,:,:,:,:,:,:,isp] = diffWY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "diffTY(" + str(species_names[isp]) + ")_mean"
    diffTY_wtsum[:,:,:,:,:,:,:,:,isp] = diffTY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "convY(" + str(species_names[isp]) + ")_mean"
    convY_wtsum[:,:,:,:,:,:,:,:,isp] = convY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
  for ir in range(0, Nr):
    fdn = "R(" + str(ir) + ")_mean"
    prr_wtsum[:,:,:,:,:,:,:,:,ir] = prr_wtsum[:,:,:,:,:,:,:,:,ir] + f["DATA"][fdn]
  print("Finished: ", fn)

#%% Temporal evolution
nt = len(fns_sorted)
hrr_t = np.zeros((nt))
rhoY_t = np.zeros((nt, Nsp))
for ifn, fn in enumerate(fns_sorted):
  f = h5py.File(fn, 'r+')
  hrr_t[ifn] = np.sum(f["DATA"]["HeatRelease_mean"])
  for isp in range(0, Nsp):
    fdn = "rhoY(" + str(species_names[isp]) + ")_mean"
    rhoY_t[ifn, isp] = np.sum(f["DATA"][fdn])
fig, ax = plt.subplots(figsize = (5, 4))
ax.plot(rhoY_t[:,N2_ID], color = "m", linestyle = "-")
ax.set_ylabel(r'$\int \textbf{HRR} \mathrm{d}V$', fontsize = labelsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='both', which='minor', labelsize=labelsize)

#%%
def get_graph(prr_bin):
  # prr_bin[1:Nr] - volume ingegral of progress rate of reaction
  # Nr - Number of reactions
  graph_res = np.zeros((Nsp, Nsp))
  for iline in range(0, transfer_matrix.shape[0]):
    ir    = int(transfer_matrix[iline, 0]) # ir-th reaction
    isp0  = int(transfer_matrix[iline, 1]) # isp0-th species
    isp1  = int(transfer_matrix[iline, 2]) # isp1-th species
    sp0   = species_names[isp0]
    sp1   = species_names[isp1]
    natom = int(transfer_matrix[iline, 3]) # Number of nitrogen atom transferred
    graph_res[isp0, isp1] = graph_res[isp0, isp1] + prr_bin[ir] * natom
  graph_res = graph_res - graph_res.T
  graph_res[graph_res<0] = 0 # Commented back
  return graph_res

#%% Test reaction rate vs. progress rate
tID = N_ID
rr = np.sum(wdot_wtsum[:,:,:,:,:,:,1,:,tID]) / (gas_CD.molecular_weights[tID]/1000)
diff = np.sum(diffY_wtsum[:,:,:,:,:,:,1,:,tID]) / (gas_CD.molecular_weights[tID]/1000)
conv = -np.sum(convY_wtsum[:,:,:,:,:,:,1,:,tID]) / (gas_CD.molecular_weights[tID]/1000)
prr_all = np.zeros((Nr))
for ir in range(0, Nr):
  prr_all[ir] = np.sum(prr_wtsum[:,:,:,:,:,:,1,:,ir],)
graph_all = get_graph(prr_all)
flux_all = -(np.sum(graph_all[tID,:]) - np.sum(graph_all[:,tID])) / gas1D.n_atoms(gas1D.species_names[tID], "N")
v3 = -(np.sum(graph_all[tID,:]) - np.sum(graph_all[:,tID])) / gas1D.n_atoms(gas1D.species_names[tID], "N")

print(r"Reaction rate [mol/s]:")
print("  By wdot(N2): ", rr)
print("  By graph   : ", v3)

print(r"Volume integral [mol/s]")
print("  diff: ", diff)
print("  conv: ", conv)
print("diff + conv: ", diff + conv)
print("flux_all", flux_all)

#%% Budget
axis_all = (0, 1, 2, 3, 4, 5, 6, 7)
axis_0123456 = (0, 1, 2, 3, 4, 5, 6)
axis_0123567 = (0, 1, 2, 3, 5, 6, 7)
axis_ZC = (0, 1, 2, 5, 6, 7)
axis_ZCF = (0, 1, 2, 5, 7)
wt_ZC_wtsum    = np.zeros((nmf, npv))
rho_ZC_wtsum    = np.zeros((nmf, npv))
rho_ZC          = np.zeros((nmf, npv))
mf_ZC_wtsum    = np.zeros((nmf, npv))
mf_ZC          = np.zeros((nmf, npv))
pv_ZC_wtsum    = np.zeros((nmf, npv))
pv_ZC          = np.zeros((nmf, npv))
rhoT_ZC_wtsum    = np.zeros((nmf, npv))
rhoT_ZC          = np.zeros((nmf, npv))
rhoY_ZC_wtsum  = np.zeros((nmf, npv, Nsp))
rhoY_ZC        = np.zeros((nmf, npv, Nsp))
wdot_ZC_wtsum  = np.zeros((nmf, npv, Nsp))
wdot_ZC        = np.zeros((nmf, npv, Nsp))
diffY_ZC_wtsum = np.zeros((nmf, npv, Nsp))
diffY_ZC       = np.zeros((nmf, npv, Nsp))
diffWY_ZC_wtsum = np.zeros((nmf, npv, Nsp))
diffWY_ZC       = np.zeros((nmf, npv, Nsp))
diffTY_ZC_wtsum = np.zeros((nmf, npv, Nsp))
diffTY_ZC       = np.zeros((nmf, npv, Nsp))
convY_ZC_wtsum = np.zeros((nmf, npv, Nsp))
convY_ZC       = np.zeros((nmf, npv, Nsp))
wt_ZC_wtsum[:] = np.sum(wt_sum, axis = axis_ZC)
rho_ZC_wtsum[:,:] = np.sum(rho_wtsum[:,:,:,:,:,:,:], axis = axis_ZC)
pv_ZC_wtsum[:,:] = np.sum(pv_wtsum[:,:,:,:,:,:,:], axis = axis_ZC)
mf_ZC_wtsum[:,:] = np.sum(mf_wtsum[:,:,:,:,:,:,:], axis = axis_ZC)
rhoT_ZC_wtsum[:,:] = np.sum(rhoT_wtsum[:,:,:,:,:,:,:], axis = axis_ZC)
for isp in range(0, Nsp):
  rhoY_ZC_wtsum[:,:,isp] = np.sum(rhoY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
  wdot_ZC_wtsum[:,:,isp] = np.sum(wdot_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
  diffY_ZC_wtsum[:,:,isp] = np.sum(diffY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
  diffWY_ZC_wtsum[:,:,isp] = np.sum(diffWY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
  diffTY_ZC_wtsum[:,:,isp] = np.sum(diffTY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
  convY_ZC_wtsum[:,:,isp] = np.sum(convY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZC)
rho_ZC[:,:] = rho_ZC_wtsum[:,:] / wt_ZC_wtsum
rhoT_ZC[:,:] = rhoT_ZC_wtsum[:,:] / wt_ZC_wtsum
mf_ZC[:,:] = mf_ZC_wtsum[:,:] / wt_ZC_wtsum
pv_ZC[:,:] = pv_ZC_wtsum[:,:] / wt_ZC_wtsum
for isp in range(0, Nsp):
  rhoY_ZC[:,:,isp] =rhoY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  wdot_ZC[:,:,isp] = wdot_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  diffY_ZC[:,:,isp] = diffY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  diffWY_ZC[:,:,isp] = diffWY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  diffTY_ZC[:,:,isp] = diffTY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
  convY_ZC[:,:,isp] = convY_ZC_wtsum[:,:,isp] / wt_ZC_wtsum
fig, ax = plt.subplots(figsize = (5, 4))
iz = 4
tID = OH_ID
print("mf 3D: ", mf_ZC[iz,:])
ax.plot(pv_ZC[iz,:], wdot_ZC[iz,:,tID], "ro-", label = r"$\rho\dot{\omega}$")
ax.plot(pv_ZC[iz,:], diffY_ZC[iz,:,tID], "bo-", label = r"$\nabla \cdot (\rho D \nabla Y)$")
ax.plot(pv_ZC[iz,:], diffWY_ZC[iz,:,tID], "ko-", label = r"$\nabla \cdot (\rho D \nabla Y)$")
ax.plot(pv_ZC[iz,:], diffTY_ZC[iz,:,tID], "mo-", label = r"$\nabla \cdot (\rho D \nabla Y)$")
ax.plot(pv_ZC[iz,:], -convY_ZC[iz,:,tID], "go-", label = r"$\nabla \cdot (\rho u Y)$")
#ax.plot(pv_ZC[iz,:], -convY_ZC[iz,:,tID]+diffY_ZC[iz,:,tID], "ko-")
ax.plot([0, 1], [0, 0])
#ax.legend(fontsize = labelsize - 8)
ymaxall = np.maximum(np.maximum(np.abs(diffY_ZC[iz,:,tID]), np.abs(-convY_ZC[iz,:,tID])), np.abs(wdot_ZC[iz,:,tID]))
ymax = np.amax(ymaxall)
ax.set_xlabel(r"$C$", fontsize = labelsize)
ax.set_ylabel(r"$\langle \textbf{Budget} | \Omega \rangle$", fontsize = labelsize)
ax.set_xlim([0, 1.05])
ax.set_ylim([-ymax, ymax])
ax.tick_params(axis='both', which='major', labelsize=24)
ax.tick_params(axis='both', which='minor', labelsize=24)
ax.set_title(r"\textbf{" + str(species_names[tID]) + "-3D}", fontsize = labelsize - 2 )
ax2 = ax.twinx()
#ax2.plot(pv_ZC[iz,:], rhoT_ZC[iz,:] / rho_ZC[iz, :], color = "m")
ax2.plot(pv_ZC[iz,:], rhoY_ZC[iz,:,tID] / rho_ZC[iz, :], color = "m")
#ax2.set_ylim(np.array([-4, 4])*1E-5)
#ax2.set_ylim([300, 2800])
#ax2.set_ylabel(r"$T~\mathrm{[K]}$", fontsize = labelsize)
ax2.set_ylabel(r"$Y_\mathrm{" + str(species_names[tID]) + "}$", fontsize = labelsize)

ax2.tick_params(axis='both', which='major', labelsize=24)
ax2.tick_params(axis='both', which='minor', labelsize=24)

#%%
prr_all = np.zeros((npv))

# All
prr_all = np.zeros((Nr))
for ir in range(0, Nr):
  prr_all[ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir],
                       axis = axis_all)
graph_all = get_graph(prr_all)
# Premixed & nonpremixed
prr_nonpre = np.zeros((Nr))
prr_premix = np.zeros((Nr))
for ir in range(0, Nr):
  prr_nonpre[ir] = np.sum(prr_wtsum[:,:,:,:,:,:,0,:,ir],
                            axis = axis_0123456)
  prr_premix[ir] = np.sum(prr_wtsum[:,:,:,4:,:,:,1,:,ir],
                            axis = axis_0123456)
graph_nonpre = get_graph(prr_nonpre)
graph_pre = get_graph(prr_premix)

#%%
# Single Z on A
prr_A = np.zeros((ntres, Nr))
for ir in range(0, Nr):
  prr_A[:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,0,:,ir], axis = (0,1,2,3,4,5))

prr_ZFA = np.zeros((nmf, nfi, ntres, Nr))
for ir in range(0, Nr):
  prr_ZFA[:,:,:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir], axis = (0, 1, 2, 4, 5))

prr_ZCFA = np.zeros((nmf, npv, nfi, ntres, Nr))
for ir in range(0, Nr):
  prr_ZCFA[:,:,:,:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir], axis = (0, 1, 2, 5))

z_A = []
path_N = []
path_N= []
path_NNH = []
path_N2O = []
imf = 4
ifi = 0
for iA in range(0, ntres):
  gt = get_graph(prr_ZFA[imf, ifi, iA, :])
  ref = np.sum(gt[N2_ID,:])
  #print("iA: ", iA, gt[N2_ID, N_ID]/ref)
  #print("iA: ", iA, gt[N2_ID, NO_ID]/ref)
  print("iA: ", iA, (gt[N2_ID, NO_ID]+gt[N2_ID, N_ID])/ref)
  print("iA: ", iA, gt[N2_ID, NNH_ID]/ref)
  print("iA: ", iA, gt[N2_ID, N2O_ID]/ref)
  print("\n")

#%% 1D test
iz1D    = 7
fn1D    = fns1D_sorted[iz1D]
df      = pd.read_csv(fn1D)
fstate  = ct.SolutionArray(gas1D)
fstate.from_pandas(df)
x_1D    = fstate.grid
T_1D    = fstate.T
rho_1D  = fstate.density
u_1D    = fstate.velocity
C_1D    = fstate.concentrations
prr_1D  = fstate.net_rates_of_progress
frr_1D  = fstate.forward_rates_of_progress
rrr_1D  = fstate.reverse_rates_of_progress
wdot_1D = fstate.net_production_rates
pv_1D = np.zeros_like(fstate.T)
for isp, spn in enumerate(gas1D.species_names):
  pv_1D = pv_1D + coeff_pv[spn] * fstate.Y[:, isp]
pv_1D = (pv_1D - pveqs0[iz1D]) / (pveqs[iz1D] - pveqs0[iz1D])
print("mf 1D: ", Zs1D[iz1D])
tID = H2_ID
Yt  = fstate.Y[:,tID]
rhoDt  = fstate.mix_diff_coeffs[:,tID] * fstate.density
grad_Yt = np.gradient(Yt, x_1D)
diff_Yt = np.gradient(rhoDt * grad_Yt, x_1D)
conv_Yt = np.gradient(rho_1D * u_1D * Yt, x_1D)
chem_Yt = wdot_1D[:,tID] * gas_CD.molecular_weights[tID]
Wbar = fstate.mean_molecular_weight / 1000
gradW_x = np.gradient(Wbar, x_1D)
phi = Yt * fstate.mix_diff_coeffs[:,tID]*fstate.density * gradW_x / Wbar
diffW_Yt = np.gradient(phi, x_1D)
DT = fstate.thermal_diff_coeffs[:,tID]
gradT_x = np.gradient(T_1D, x_1D)
diffT_Yt = DT * gradT_x / T_1D
diffT_Yt = np.gradient(diffT_Yt, x_1D)
fig, ax = plt.subplots(figsize = (5, 4))
ax.plot(pv_1D, chem_Yt, "ro-", label = r"$\rho\dot{\omega}$")
ax.plot(pv_1D, diff_Yt, "bo-", label = r"$\nabla \cdot (\rho D \nabla Y)$")
ax.plot(pv_1D, diffW_Yt, "ko-", label = r"$\nabla \cdot (\rho D Y \frac{\nabla \overline{W}}{\overline{W}})$")
ax.plot(pv_1D, diffT_Yt, "mo-", label = r"$\nabla \cdot (\frac{\rho D_T}{T} \nabla T)$")
ax.plot(pv_1D, -conv_Yt, "go-", label = r"$\nabla \cdot (\rho u Y)$")

ax.plot([0, 1], [0, 0])
ax.legend(fontsize = labelsize - 8,
          bbox_to_anchor=(0.8, 1.5))
ymaxall = np.maximum(np.maximum(np.abs(diffY_ZC[iz,:,tID]), np.abs(-convY_ZC[iz,:,tID])), np.abs(wdot_ZC[iz,:,tID]))
ymax = np.amax(ymaxall)
ax.set_xlabel(r"$C$", fontsize = labelsize)
ax.set_ylabel(r"$ \textbf{Budget} $", fontsize = labelsize)
ax.set_xlim([0.0, 1.05])
#ax.set_ylim([-ymax, ymax])
ax.tick_params(axis='both', which='major', labelsize=24)
ax.tick_params(axis='both', which='minor', labelsize=24)
ax.set_title(r"\textbf{" + str(species_names[tID]) + "-1D}", fontsize = labelsize - 2 )
ax2 = ax.twinx()
ax2.plot(pv_1D, Yt, color = "m")
#ax2.set_ylim(np.array([-4, 4])*1E-5)
ax2.set_ylabel(r"$Y_\mathrm{" + str(species_names[tID]) + "}$", fontsize = labelsize)
ax2.tick_params(axis='both', which='major', labelsize=24)
ax2.tick_params(axis='both', which='minor', labelsize=24)

#fig, ax = plt.subplots(figsize = (5, 4))
#ax.plot(pv_1D, fstate.fstate.thermal_diff_coeffs)

#%%
tID     = N2_ID
dx      = np.zeros_like(x_1D)
dx[0:-1] = x_1D[1:] - x_1D[0:-1]
dx[-1]  = dx[-2]
# 1
v1 = u_1D[-1]*C_1D[-1,tID] - u_1D[0]*C_1D[0,tID]
# 2
wdot_int = np.sum(wdot_1D[:,tID] * dx)
v2 = wdot_int
# 3
prr_wtsum_1D  = np.zeros((Nr))
for ir in range(0, Nr):
  prr_wtsum_1D[ir] = np.sum(dx * prr_1D[:,ir])
graph_1D = get_graph(prr_wtsum_1D)
flux1D = np.sum(gt[N2_ID,:])
print("Reference flux 1D", flux1D)

#from_others = np.sum(graph_1D[:,tID]) - graph_1D[tID,tID]
#to_others = np.sum(graph_1D[tID,:])
#v3 = (to_others - from_others)/2.
v3 = -(np.sum(graph_1D[tID,:]) - np.sum(graph_1D[:,tID])) / gas1D.n_atoms(gas1D.species_names[tID], "N")

print("1. Outflow flux:", v1, "kmol/m2/s")
print("2. Production of NO:", v2, "kmol/m2/s")
print("3. Production of NO from transfer matrix:", v3, "kmol/m2/s")
fig, ax = plt.subplots(figsize = (5, 4))
ax.plot(fstate.grid, fstate.Y[:,NO_ID]/np.amax(fstate.Y[:,NO_ID]), color="b", linestyle = "-", label = "Y(NO)")
ax.plot(fstate.grid, fstate.T/np.amax(fstate.T), color="r", linestyle = "-", label = "T")
ax.set_xlim([0.029, 0.032])

import graphviz
fg_1D = graphviz.Digraph(filename="part_premix_rich",
                           format="png",
                          graph_attr={"nodesep" : "0.0",
                                      "pad" : "0.1",
                                      "fontname" : "times bold"},
)
gt = graph_1D
nodesize = 50
labelsize = 40
flux_ref = np.sum(gt[N2_ID,:])

spg_part = ["N2", "NO", "N", "NNH", "N2O"]
xref = np.array([0, 1, 2, 3]);
yref = np.array([0, 1, 2, 3])
pos = {}
fg_1D.attr(size=f"{len(spg_part)}")
for i in range(0, len(spg_part)):
  nnode = f"{spg_part[i]}"
  fg_1D.node(nnode, label = nnode, fontsize=f"{nodesize}")

for i in range(0, len(spg_part)):
  for j in range(0, len(spg_part)):
    spo = spg_part[i]
    isp0 = gas1D.species_index(spo)
    origin = f"{spo}"

    spd = spg_part[j]
    isp1 = gas1D.species_index(spd)
    destin = f"{spd}"

    color='black'
    if ((spo == "N2" and spd == "N") or
        spo == "N" and spd == "N2"):
      color = 'red'
    if ((spo == "N2" and spd == "NO") or
        spo == "NO" and spd == "N2"):
      color = 'red'
    if ((spo == "N2" and spd == "NNH") or
        spo == "NNH" and spd == "N2"):
      color = 'purple'
    if ((spo == "N2" and spd == "N2O") or
        spo == "N2O" and spd == "N2"):
      color = 'blue'


    width = 40*gt[isp0, isp1]/flux_ref
    width = max(1, width)

    label = f"  {100*gt[isp0, isp1]/flux_ref:.2f} %\n\n"
    add_edge = (gt[isp0, isp1] / flux_ref) > 0.01
    if add_edge:
      fg_1D.edge(origin, destin,
                   label=label,
                   penwidth=str(width),
                   color=color,
                   fontcolor=color,
                   arrowhead="vee",
                   fontsize = f"{labelsize}")
fg_1D.render()
fg_1D

#%%
import graphviz
fg_part = graphviz.Digraph(filename="part_premix_rich",
                           format="png",
                          graph_attr={"nodesep" : "0.0",
                                      "pad" : "0.1",
                                      "fontname" : "times bold"},
)
gt = graph_nonpre
nodesize = 50
labelsize = 40
flux_ref = np.sum(gt[N2_ID,:])
flux_ref = np.sum(gt[N2_ID,:]) - np.sum(gt[:,N2_ID])

spg_part = ["N2", "NO", "N", "NNH", "N2O", "NH", "NH2", "NO2",]
xref = np.array([0, 1, 2, 3]);
yref = np.array([0, 1, 2, 3])
pos = {}
fg_part.attr(size=f"{len(spg_part)}")
for i in range(0, len(spg_part)):
  nnode = f"{spg_part[i]}"
  fg_part.node(nnode, label = nnode, fontsize=f"{nodesize}")

for i in range(0, len(spg_part)):
  for j in range(0, len(spg_part)):
    spo = spg_part[i]
    isp0 = gas1D.species_index(spo)
    origin = f"{spo}"

    spd = spg_part[j]
    isp1 = gas1D.species_index(spd)
    destin = f"{spd}"

    color='black'
    if ((spo == "N2" and spd == "N") or
        spo == "N" and spd == "N2"):
      color = 'red'
    if ((spo == "N2" and spd == "NO") or
        spo == "NO" and spd == "N2"):
      color = 'red'
    if ((spo == "N2" and spd == "NNH") or
        spo == "NNH" and spd == "N2"):
      color = 'purple'
    if ((spo == "N2" and spd == "N2O") or
        spo == "N2O" and spd == "N2"):
      color = 'blue'


    width = 40*gt[isp0, isp1]/flux_ref
    width = max(1, width)

    label = f"  {100*gt[isp0, isp1]/flux_ref:.2f} %\n\n"
    add_edge = (gt[isp0, isp1] / flux_ref) > 0.002
    if add_edge:
      fg_part.edge(origin, destin,
                   label=label,
                   penwidth=str(width),
                   color=color,
                   fontcolor=color,
                   arrowhead="vee",
                   fontsize = f"{labelsize}")
fg_part.render()
fg_part
#%%
axis_ZCF = (0, 1, 2, 5, 7)
wt_ZCF_wtsum    = np.zeros((nmf, npv, nfi))
wdot_ZCF_wtsum  = np.zeros((nmf, npv, nfi, Nsp))
wdot_ZCF        = np.zeros((nmf, npv, nfi, Nsp))
diffY_ZCF_wtsum = np.zeros((nmf, npv, nfi, Nsp))
diffY_ZCF       = np.zeros((nmf, npv, nfi, Nsp))
convY_ZCF_wtsum = np.zeros((nmf, npv, nfi, Nsp))
convY_ZCF       = np.zeros((nmf, npv, nfi, Nsp))

wt_ZCF_wtsum[:] = np.sum(wt_sum, axis = axis_ZCF)
for isp in range(0, Nsp):
  wdot_ZCF_wtsum[:,:,:,isp] = np.sum(wdot_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZCF)
  diffY_ZCF_wtsum[:,:,:,isp] = np.sum(diffY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZCF)
  convY_ZCF_wtsum[:,:,:,isp] = np.sum(convY_wtsum[:,:,:,:,:,:,:,:,isp], axis = axis_ZCF)

for isp in range(0, Nsp):
  wdot_ZCF[:,:,:,isp] = wdot_ZCF_wtsum[:,:,:,isp] / wt_ZCF_wtsum
  diffY_ZCF[:,:,:,isp] = diffY_ZCF_wtsum[:,:,:,isp] / wt_ZCF_wtsum
  convY_ZCF[:,:,:,isp] = convY_ZCF_wtsum[:,:,:,isp] / wt_ZCF_wtsum

#%%
tID = N2O_ID
sum_chem = np.sum(wdot_ZCF_wtsum[4:,:,1,tID]) / (gas_CD.molecular_weights[tID]/1000) / gas1D.n_atoms(gas1D.species_names[tID], "N")
print("react: ", sum_chem)
sum_diff = np.sum(diffY_ZCF_wtsum[4:,:,1,tID]) / (gas_CD.molecular_weights[tID]/1000) / gas1D.n_atoms(gas1D.species_names[tID], "N")
print("diff: ", sum_diff)
sum_conv = np.sum(convY_ZCF_wtsum[4:,:,1,tID]) / (gas_CD.molecular_weights[tID]/1000) / gas1D.n_atoms(gas1D.species_names[tID], "N")
print("conv: ", -sum_conv)
print("diff+conv: ", sum_diff-sum_conv)
print("flux_ref: ", flux_ref)
print("tran/flux_ref", (sum_diff-sum_conv)/flux_ref*100, "%")

#for isp in range(0, gas_CD.n_species):
#  print(gas_CD.species_name(tID)+"->"+gas_CD.species_name(isp)+":", gt[tID, isp]/flux_ref)

#%% Dependence of RPA on post-flame residence time
prr_A = np.zeros((ntres, Nr))
for ir in range(0, Nr):
  prr_A[:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,0,:,ir], axis = (0,1,2,3,4,5))

prr_ZFA = np.zeros((nmf, nfi, ntres, Nr))
for ir in range(0, Nr):
  prr_ZFA[:,:,:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir], axis = (0, 1, 2, 4, 5))

prr_ZCFA = np.zeros((nmf, npv, nfi, ntres, Nr))
for ir in range(0, Nr):
  prr_ZCFA[:,:,:,:,ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir], axis = (0, 1, 2, 5))


lw = 2.5; markersize = 12; markeredgewidth = 3;
fontsize = 20
z_A = []
path_N = []
path_N= []
path_NNH = []
path_N2O = []
# (10, 19)
imf = 10; smf = slice(10, 19)
ifi = 0
zmin = np.linspace(0, 0.1, 20)[smf.start]
zmax = np.linspace(0, 0.1, 20)[smf.stop]
zst = 0.0252
strt = r"$Z /Z_\mathrm{st}\in [" + "%.1f"%(zmin/zst) + ",~"
strt = strt + "%.1f"%(zmax/zst) + "] $"

res_A = np.linspace(0.05E-4, 1.2E-4, ntres)
res_nonpre = np.zeros((ntres, 3))
res_pre = np.zeros((ntres, 3))
for iA in range(0, ntres):
  prr_local = np.sum(prr_ZFA[smf, 0, iA, :], axis = (0))
  gt = get_graph(prr_local)
  ref = np.sum(gt[N2_ID,:])
  r0 = (gt[N2_ID, NO_ID]+gt[N2_ID, N_ID])/ref
  res_nonpre[iA, 0] = r0
  r1 = gt[N2_ID, NNH_ID]/ref
  res_nonpre[iA, 1] = r1
  r2 = gt[N2_ID, N2O_ID]/ref
  res_nonpre[iA, 2] = r2

  prr_local = np.sum(prr_ZFA[smf, 1, iA, :], axis = (0))
  gt = get_graph(prr_local)
  ref = np.sum(gt[N2_ID,:])
  r0 = (gt[N2_ID, NO_ID]+gt[N2_ID, N_ID])/ref
  res_pre[iA, 0] = r0
  r1 = gt[N2_ID, NNH_ID]/ref
  res_pre[iA, 1] = r1
  r2 = gt[N2_ID, N2O_ID]/ref
  res_pre[iA, 2] = r2


res_nonpre = pd.DataFrame({ 'A' : res_A[:],
                      'thermal': res_nonpre[:,0],
                      'NNH': res_nonpre[:,1],
                      'N2O': res_nonpre[:,2]})
res_pre = pd.DataFrame({ 'A' : res_A[:],
                      'thermal': res_pre[:,0],
                      'NNH': res_pre[:,1],
                      'N2O': res_pre[:,2]})

res_nonpre.dropna(subset=['thermal'])
res_pre.dropna(subset=['thermal'])

A_nonpre = res_nonpre['A'].values * 10000.
nonpre0 = res_nonpre['thermal'].values
nonpre1 = res_nonpre['NNH'].values
nonpre2 = res_nonpre['N2O'].values
non_pre = res_nonpre['A'].values
A_pre = res_pre['A'].values * 10000.
pre0 = res_pre['thermal'].values
pre1 = res_pre['NNH'].values
pre2 = res_pre['N2O'].values

figunit_x = 3.0; figunit_y = 3.0
fig, axs = plt.subplots(figsize = (figunit_x, figunit_y*2), ncols = 1, nrows = 2)
for iplot in range(0, 2):
  ax = axs[iplot]
  if (iplot == 0):
    ax.plot(A_nonpre, nonpre0, color = "r", linewidth = 2.5,
        marker = "o", markersize = markersize,
        markerfacecolor = "w",
        markeredgecolor = "r", markeredgewidth = markeredgewidth)
    ax.set_title(strt, fontsize = 24)
    ax.plot(A_nonpre, nonpre1, color = "purple", linewidth = 2.5,
        marker = "o", markersize = markersize,
        markerfacecolor = "w",
        markeredgecolor = "purple", markeredgewidth = markeredgewidth)
    ax.plot(A_nonpre, nonpre2, color = "blue", linewidth = 2.5,
        marker = "o", markersize = markersize,
        markerfacecolor = "w",
        markeredgecolor = "blue", markeredgewidth = markeredgewidth)
  elif (iplot == 1):
    ax.plot(A_pre, pre0, color = "r", linewidth = 2.5,
        marker = "s", markersize = markersize,
        markerfacecolor = "w",
        markeredgecolor = "r", markeredgewidth = markeredgewidth)
    ax.plot(A_pre, pre1, color = "purple", linewidth = 2.5,
        marker = "s", markersize = markersize, alpha = 1.0,
        markerfacecolor = "w",
        markeredgecolor = "purple", markeredgewidth = markeredgewidth)
    ax.plot(A_pre, pre2, color = "blue", linewidth = 2.5,
        marker = "s", markersize = markersize,
        markerfacecolor = "w",
        markeredgecolor = "blue", markeredgewidth = markeredgewidth)
  if iplot == 1:
    ax.set_xlabel(r'$\alpha_\mathrm{p} \times \mathrm{10^{-4}~[s]}$', fontsize = 20)
  #ax.set_ylabel("", fontsize = 24)
  if iplot == 1:
    ax.set_xticks([0.6, 1.2])
  else:
    ax.set_xticks([0.6, 1.2])
    ax.set_xticklabels([])
  ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
  ax.set_yticklabels([r"$0.0$", " ", r"$0.5$", " ", r"$1.0$"])
  ax.tick_params(axis='both', which='major', labelsize=24)
  ax.tick_params(axis='both', which='minor', labelsize=24)
  ax.grid()
pltn = "./RPAA_Z=" + str(smf.start) + "-" + str(smf.stop) + ".png"
plt.savefig(pltn, dpi = 400, bbox_inches = "tight")

#%%