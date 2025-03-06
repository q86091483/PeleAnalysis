#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
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

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
labelsize = 16

# Input
# Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_NOx/"

# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_14*.h5")

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
ntres = tet.shape[7]
zout = np.linspace(2.5E-3, 2.575E-1, nmf)
pvout = np.linspace(300, 2800, npv)

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
idx_H2  = gas1D.species_index("H2")
idx_O2  = gas1D.species_index("O2")
idx_H2O = gas1D.species_index("H2O")
idx_OH  = gas1D.species_index("OH")
idx_H   = gas1D.species_index("H")
idx_O   = gas1D.species_index("O")
idx_N   = gas1D.species_index("N")
idx_N2  = gas1D.species_index("N2")
idx_NO  = gas1D.species_index("NO")
idx_NO2 = gas1D.species_index("NO2")
idx_N2O = gas1D.species_index("N2O")
idx_NNH = gas1D.species_index("NNH")

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

# x, y, z, Z, C, zone, FI, tres
dim_cond = (nx, ny, nz, nmf, npv, nzo, nfi, ntres)
dim_cond_sp = dim_cond + (Nsp,)
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

for ifn, fn in enumerate(fns):
  f = h5py.File(fn, 'r+')
  wt_sum  = wt_sum        + f["DATA"]["volume_mean"]
  rho_wtsum = rho_wtsum   + f["DATA"]["rho_mean"]
  rhoT_wtsum = rhoT_wtsum + f["DATA"]["rhoT_mean"]
  rhoT2_wtsum = rhoT2_wtsum + f["DATA"]["rhoT2_mean"]
  agepv_wtsum = agepv_wtsum + f["DATA"]["agepv_1_mean"]
  mf_wtsum = mf_wtsum     + f["DATA"]["mixture_fraction_mean"]
  pv_wtsum = pv_wtsum     + f["DATA"]["pv_mean"]
  hrr_wtsum = hrr_wtsum   + f["DATA"]["HeatRelease_mean"]
  hrr2_wtsum = hrr2_wtsum   + f["DATA"]["HeatRelease2_mean"]
  for isp in range(0, Nsp):
    fdn = "rhoY(" + str(species_names[isp]) + ")_mean"
    rhoY_wtsum[:,:,:,:,:,:,:,:,isp] = rhoY_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "rhoY2(" + str(species_names[isp]) + ")_mean"
    rhoY2_wtsum[:,:,:,:,:,:,:,:,isp] = rhoY2_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "wdot(" + str(species_names[isp]) + ")_mean"
    wdot_wtsum[:,:,:,:,:,:,:,:,isp] = wdot_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
    fdn = "wdot2(" + str(species_names[isp]) + ")_mean"
    wdot2_wtsum[:,:,:,:,:,:,:,:,isp] = wdot2_wtsum[:,:,:,:,:,:,:,:,isp] + f["DATA"][fdn]
  print("Finished: ", fn)
#%%
# 0 - x
# 1 - y
# 2 - z
# 3 - Z
# 4 - C
# 5 - zone
# 6 - FI
# 7 - tres
# <|Z,C,A>
axis_sum = (0, 1, 2, 5, 6)
wt_ZCA_wtsum     = np.sum(wt_sum, axis=axis_sum)
rho_ZCA_wtsum    = np.sum(rho_wtsum, axis=axis_sum)
mf_ZCA_wtsum     = np.sum(mf_wtsum, axis=axis_sum)
pv_ZCA_wtsum     = np.sum(pv_wtsum, axis=axis_sum)
rhoT_ZCA_wtsum   = np.sum(rhoT_wtsum, axis=axis_sum)
rhoT2_ZCA_wtsum  = np.sum(rhoT2_wtsum, axis=axis_sum)
agepv_ZCA_wtsum   = np.sum(agepv_wtsum, axis=axis_sum)
hrr_ZCA_wtsum    = np.sum(hrr_wtsum, axis=axis_sum)
hrr2_ZCA_wtsum   = np.sum(hrr2_wtsum, axis=axis_sum)
rhoY_ZCA_wtsum   = np.sum(rhoY_wtsum, axis=axis_sum)
rhoY2_ZCA_wtsum  = np.sum(rhoY2_wtsum, axis=axis_sum)
wdot_ZCA_wtsum   = np.sum(wdot_wtsum, axis=axis_sum)
wdot2_ZCA_wtsum  = np.sum(wdot2_wtsum, axis=axis_sum)

rho_ZCA    = rho_ZCA_wtsum / wt_ZCA_wtsum
mf_ZCA     = mf_ZCA_wtsum / wt_ZCA_wtsum
pv_ZCA     = pv_ZCA_wtsum / wt_ZCA_wtsum
rhoT_ZCA   = rhoT_ZCA_wtsum / wt_ZCA_wtsum
rhoT2_ZCA  = rhoT2_ZCA_wtsum / wt_ZCA_wtsum
agepv_ZCA  = agepv_ZCA_wtsum / wt_ZCA_wtsum
hrr_ZCA    = hrr_ZCA_wtsum / wt_ZCA_wtsum
hrr2_ZCA   = hrr2_ZCA_wtsum / wt_ZCA_wtsum
rhoY_ZCA   = np.zeros_like(rhoY_ZCA_wtsum)
rhoY2_ZCA  = np.zeros_like(rhoY2_ZCA_wtsum)
wdot_ZCA   = np.zeros_like(wdot_ZCA_wtsum)
wdot2_ZCA  = np.zeros_like(wdot2_ZCA_wtsum)
for isp in range(0, Nsp):
  rhoY_ZCA[:,:,:,isp] = rhoY_ZCA_wtsum[:,:,:,isp] / wt_ZCA_wtsum
  rhoY2_ZCA[:,:,:,isp] = rhoY2_ZCA_wtsum[:,:,:,isp] / wt_ZCA_wtsum
  wdot_ZCA[:,:,:,isp] = wdot_ZCA_wtsum[:,:,:,isp] / wt_ZCA_wtsum
  wdot2_ZCA[:,:,:,isp] = wdot2_ZCA_wtsum[:,:,:,isp] / wt_ZCA_wtsum

#%% Analysis
# Single Z for integrating <Y(NO)|Z,C,A>
Zt      = 0.05
iz1D    = np.argmax((Zs1D > Zt))-1 ; print("In 1D: Z0 = ", Zs1D[iz1D])
iz1D_plus = iz1D + 1; print("In 1D: Z1 = ", Zs1D[iz1D_plus])
iz3D    = np.argmax((zout > Zt)) ; print("In 3D: Z = ", zout[iz3D])
# 1D - C, T, Y
fn1D      = fns1D_sorted[iz1D]
fn1D_plus = fns1D_sorted[iz1D_plus]

df        = pd.read_csv(fn1D)
df_plus   = pd.read_csv(fn1D_plus)

fstate  = ct.SolutionArray(gas1D)
fstate.from_pandas(df)
fstate_plus  = ct.SolutionArray(gas1D)
fstate_plus.from_pandas(df_plus)

net_production_rates      = fstate.net_production_rates
net_production_rates_plus = fstate_plus.net_production_rates

pv1D        = np.zeros_like(fstate.T)
pv1D_plus   = np.zeros_like(fstate_plus.T)
for isp, spn in enumerate(gas1D.species_names):
  pv1D      = pv1D + coeff_pv[spn] * fstate.Y[:, isp]
  pv1D_plus = pv1D_plus + coeff_pv[spn] * fstate_plus.Y[:, isp]

pv1D        = (pv1D - pveqs0[iz1D]) / (pveqs[iz1D] - pveqs0[iz1D])
pv1D_plus   = (pv1D_plus - pveqs0[iz1D]) / (pveqs[iz1D] - pveqs0[iz1D])

T_1D        = fstate.T[:]
T_1D_plus   = fstate_plus.T[:]

Y_1D        = np.zeros((T_1D.shape[0], Nsp))
Y_1D_plus   = np.zeros((T_1D_plus.shape[0], Nsp))

for isp in range(0, Nsp):
  Y_1D[:, isp]      = fstate.Y[: ,isp]
  Y_1D_plus[:, isp] = fstate_plus.Y[: ,isp]

# 1D - residence time
indx_pf     = (fstate.T[:] > 1750.)
x_pf_1D     = fstate.grid[indx_pf]
u_pf_1D     = fstate.velocity[indx_pf]
agepv_1D    = np.zeros_like(u_pf_1D)
for i in range(1, agepv_1D.shape[0]):
  agepv_1D[i] = agepv_1D[i-1] + (x_pf_1D[i] - x_pf_1D[i-1]) / u_pf_1D[i]

indx_pf_plus     = (fstate_plus.T[:] > 1750.)
x_pf_1D_plus     = fstate_plus.grid[indx_pf_plus]
u_pf_1D_plus     = fstate_plus.velocity[indx_pf_plus]
agepv_1D_plus    = np.zeros_like(u_pf_1D_plus)
for i in range(1, agepv_1D_plus.shape[0]):
  agepv_1D_plus[i] = agepv_1D_plus[i-1] + (x_pf_1D_plus[i] - x_pf_1D_plus[i-1]) / u_pf_1D_plus[i]
# 3D - <Y(NO)|A> at single Z
axis_C = (0)
wt_A_wtsum    = np.sum(wt_ZCA_wtsum[iz3D,:,:], axis = axis_C)
rho_A_wtsum   = np.sum(rho_ZCA_wtsum[iz3D,:,:], axis = axis_C)
pv_A_wtsum    = np.sum(pv_ZCA_wtsum[iz3D,:,:], axis = axis_C)
rhoT_A_wtsum  = np.sum(rhoT_ZCA_wtsum[iz3D,:,:], axis = axis_C)
agepv_A_wtsum = np.sum(agepv_ZCA_wtsum[iz3D,:,:], axis = axis_C)
sl = (wt_A_wtsum > 0)
rho_A      = rho_A_wtsum[sl] / wt_A_wtsum[sl]
rhoT_A     = rhoT_A_wtsum[sl] / wt_A_wtsum[sl]
agepv_A    = agepv_A_wtsum[sl] / wt_A_wtsum[sl]
pv_A       = pv_A_wtsum[sl] / wt_A_wtsum[sl]
wdot_A     = np.zeros((rho_A.shape[0], Nsp))
rhoY_A     = np.zeros((rho_A.shape[0], Nsp))
rhoY2_A     = np.zeros((rho_A.shape[0], Nsp))
Y_A        = np.zeros_like(rhoY_A)
for isp in range(0, Nsp):
  wdot_A[:, isp] = np.sum(wdot_ZCA_wtsum[iz3D,:,:,isp], axis = axis_C)[sl] / wt_A_wtsum[sl]
  rhoY_A[:, isp] = np.sum(rhoY_ZCA_wtsum[iz3D,:,:,isp], axis = axis_C)[sl] / wt_A_wtsum[sl]
  rhoY2_A[:, isp] = np.sum(rhoY2_ZCA_wtsum[iz3D,:,:,isp], axis = axis_C)[sl] / wt_A_wtsum[sl]
  Y_A[:, isp] = rhoY_A[:, isp] / rho_A[:]
T_A = rhoT_A / rho_A

# 3D - <Y(NO)|C> at single Z
axis_A = (1)
wt_C_wtsum    = np.sum(wt_ZCA_wtsum[iz3D,:,:], axis = axis_A)
rho_C_wtsum   = np.sum(rho_ZCA_wtsum[iz3D,:,:], axis = axis_A)
pv_C_wtsum    = np.sum(pv_ZCA_wtsum[iz3D,:,:], axis = axis_A)
rhoT_C_wtsum  = np.sum(rhoT_ZCA_wtsum[iz3D,:,:], axis = axis_A)
agepv_C_wtsum  = np.sum(agepv_ZCA_wtsum[iz3D,:,:], axis = axis_A)
sl = (wt_C_wtsum > 0)
rho_C      = rho_C_wtsum[sl] / wt_C_wtsum[sl]
rhoT_C     = rhoT_C_wtsum[sl] / wt_C_wtsum[sl]
agepv_C    = agepv_C_wtsum[sl] / wt_C_wtsum[sl]
pv_C       = pv_C_wtsum[sl] / wt_C_wtsum[sl]
wdot_C     = np.zeros((rho_C.shape[0], Nsp))
rhoY_C     = np.zeros((rho_C.shape[0], Nsp))
rhoY2_C    = np.zeros((rho_C.shape[0], Nsp))
Y_C        = np.zeros_like(rhoY_C)
for isp in range(0, Nsp):
  rhoY_C[:, isp] = np.sum(rhoY_ZCA_wtsum[iz3D,:,:,isp], axis = axis_A)[sl] / wt_C_wtsum[sl]
  rhoY2_C[:, isp] = np.sum(rhoY2_ZCA_wtsum[iz3D,:,:,isp], axis = axis_A)[sl] / wt_C_wtsum[sl]
  wdot_C[:, isp] = np.sum(wdot_ZCA_wtsum[iz3D,:,:,isp], axis = axis_A)[sl] / wt_C_wtsum[sl]
  Y_C[:, isp] = rhoY_C[:, isp] / rho_C
T_C = rhoT_C / rho_C

#%%
figunit_x = 2; npx = 3
figunit_y = 2; npy = 1
fig, axs = plt.subplots(figsize = (npx * figunit_x+3, npy * figunit_y),
                        ncols = npx, nrows = npy)
lw = 3.5
for iplot in range(0, npx):
  ax = axs[iplot]
  c0 = "teal"; lst0 = "--"
  c1 = "crimson"; lst1 = "--"
  if iplot == 0:
    ax.plot(pv1D, Y_1D[:, idx_NO]*1E4, linestyle=lst0, color=c0, linewidth = lw)
    ax.plot(pv1D_plus, Y_1D_plus[:, idx_NO]*1E4, linestyle=lst1, color=c1, linewidth = lw)
    ax.plot(pv_C, Y_C[:, idx_NO]*1E4, linestyle="-", color="k", linewidth = lw)
    rhoY2_fms_C = rhoY2_C[:,idx_NO] - rho_C[:]*Y_C[:, idx_NO]*Y_C[:, idx_NO]
    Y_fms = np.sqrt(rhoY2_fms_C / rho_C)
    y0 = Y_C[:, idx_NO] - Y_fms
    y1 = Y_C[:, idx_NO] + Y_fms
    ax.fill_between(pv_C, y0*1E4, y1*1E4, color="gray", alpha=0.5)
    ax.set_xlim([0, 1.2])
    ax.set_ylim([0, 1.5])
    ax.set_xlabel(r"$C$", fontsize = 24)
    ax.set_title(r'$Y_\mathrm{NO} \times 10^{-4}$', fontsize = 24)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(3,1))

    #ax.plot([0.05, 0.3], [1.1, 1.1], linestyle="-", color="k", linewidth = lw-1.0)
    #ax.fill_between([0.0, 0.3], [1.2, 1.2], [1.0,1.0], color="gray", alpha=0.5)
    ax.text(0.05, 1.0, r'$Z_0=0.05$', color = "black", fontsize = 28)

  if iplot == 1:
    ax.plot(pv1D[indx_pf], agepv_1D*1E4, linestyle=lst0, color=c0, linewidth = lw)
    ax.plot(pv1D_plus[indx_pf_plus], agepv_1D_plus*1E4, linestyle=lst1, color=c1, linewidth = lw)
    ax.plot(pv_C, agepv_C[:]*1E4, linestyle="-", color="k", linewidth = lw)
    ax.set_xlim([0.8, 1.06]); ax.set_xticks([0.9, 1.0])
    ax.set_ylim([0, 1]); ax.set_yticks([0, 1])
    ax.set_xlabel(r"$C$", fontsize = 24)
    ax.set_title(r'$\alpha_\mathrm{p} \times 10^{-4}~\mathrm{[s]}$', fontsize = 24)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(3,1))
  elif iplot == 2:
    ax.plot(agepv_1D*1E4, fstate.Y[indx_pf,idx_NO]*1E4, linestyle=lst0, color=c0, linewidth = lw)
    ax.plot(agepv_1D_plus*1E4, fstate_plus.Y[indx_pf_plus,idx_NO]*1E4, linestyle=lst1, color=c1, linewidth = lw)
    ax.plot(agepv_A*1E4, Y_A[:, idx_NO]*1E4, linestyle="-", color="k", linewidth = lw)
    rhoY2_fms_A = rhoY2_A[:,idx_NO] - rho_A[:]*Y_A[:, idx_NO]*Y_A[:, idx_NO]
    Y_fms = np.sqrt(rhoY2_fms_A / rho_A)
    y0 = Y_A[:, idx_NO] - Y_fms
    y1 = Y_A[:, idx_NO] + Y_fms
    ax.fill_between(agepv_A*1E4, y0*1E4, y1*1E4, color="gray", alpha=0.5)
    ax.set_xlim([0, 1.2]); ax.set_xticks([0, 0.6, 1.2])
    ax.set_ylim([0, 1.5])
    ax.set_xlabel(r'$\alpha_\mathrm{p} \times 10^{-4} ~\mathrm{[s]}$', fontsize = 24)
    ax.set_title(r'$Y_\mathrm{NO} \times 10^{-4}$', fontsize = 24)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(3,1))

  ax.tick_params(axis='both', which='major', labelsize=24)
  ax.tick_params(axis='both', which='minor', labelsize=24)
plt.subplots_adjust(hspace = 0.3)
plt.savefig("./YNOA-" + str(int(Zt*100))+".png", dpi = 400, bbox_inches = "tight")
#%%

