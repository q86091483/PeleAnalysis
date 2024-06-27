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
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix_FI/"
# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_1*_Y.h5")
# Data folder - where 1D flame data is stored
data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/UnstrainedPremixed1D"
# Derived parameter
fn = fns[0]
f = h5py.File(fn, 'r+')
tet = np.array(f["DATA"]["volume_mean"])
nx = tet.shape[0]
ny = tet.shape[1]
nz = tet.shape[2]
nmf = tet.shape[3]
npv = tet.shape[4]

zout = np.linspace(0, 1, nmf)
pvout = np.linspace(0, 1, npv)

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
Zs = []
for ifn, fn in enumerate(fns1D_sorted):
  Zs.append(get_Z_key(fn))

# Domain size
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

# Load 0D states
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
# Load 1D states
fn_CD = "/scratch/b/bsavard/zisen347/PeleAnalysis/RJICF/Data/CounterDiffusion1D/a=2.61E+03.csv"
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
H2_ID = gas_CD.species_index("H2")
#%%
# xyzZC
wt_xyzZC_sum = np.zeros((nx, ny, nz, nmf, npv))
rho_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))
T_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))
mf_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))
pv_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))
hrr_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))
rhoY_H2_xyzZC_wtsum = np.zeros((nx, ny, nz, nmf, npv))


for ifn, fn in enumerate(fns):
  f = h5py.File(fn, 'r+')
  print(fn)
  wt_xyzZC_sum  = wt_xyzZC_sum  + f["DATA"]["volume_mean"]
  rho_xyzZC_wtsum = rho_xyzZC_wtsum + f["DATA"]["rho_mean"]
  T_xyzZC_wtsum = T_xyzZC_wtsum + f["DATA"]["temp_mean"]
  mf_xyzZC_wtsum = mf_xyzZC_wtsum + f["DATA"]["mixture_fraction_mean"]
  pv_xyzZC_wtsum = pv_xyzZC_wtsum + f["DATA"]["pv_mean"]
  hrr_xyzZC_wtsum = hrr_xyzZC_wtsum + f["DATA"]["HeatRelease_mean"]
  rhoY_H2_xyzZC_wtsum = rhoY_H2_xyzZC_wtsum + f["DATA"]["rhoY(H2)_mean"]
#%%

axis_sum = (0, 1, 2, 4)
wt_Z_wtsum = np.sum(wt_xyzZC_sum, axis=axis_sum)
rho_Z_wtsum = np.sum(rho_xyzZC_wtsum, axis=axis_sum)
T_Z_wtsum = np.sum(T_xyzZC_wtsum, axis=axis_sum)
mf_Z_wtsum = np.sum(mf_xyzZC_wtsum, axis=axis_sum)
hrr_Z_wtsum = np.sum(hrr_xyzZC_wtsum, axis=axis_sum)
rhoY_H2_Z_wtsum = np.sum(rhoY_H2_xyzZC_wtsum, axis=axis_sum)


#%%
mf_Z = mf_Z_wtsum / wt_Z_wtsum 
T_Z = T_Z_wtsum / wt_Z_wtsum 
rho_Z = rho_Z_wtsum / wt_Z_wtsum 
rhoY_H2_Z = rhoY_H2_Z_wtsum / wt_Z_wtsum

fig, ax = plt.subplots()
ax.plot(mf_CD, fstate_CD.T, label=r"$1D-90\% \chi_\mathrm{st}^{extinction}$") 
ax.plot(zout, T_Z, label=r"$\langle T | Z \rangle-2D$")
ax.plot(zout, states_0D_min.T, color="k", linestyle="--", label="Pure mix")
ax.plot(zout, states_0D_max.T, color="k", linestyle="-", label="Equilibirum")
ax.legend()
ax.legend(fontsize=20)
#ax.set_ylim([300, 2500])
ax.set_ylabel(r"$\langle Z|Z \rangle$", fontsize=20)

fig, ax = plt.subplots()
ax.plot(mf_CD, fstate_CD.Y[:,H2_ID], color="b", linewidth=3.0, linestyle = "--",
        label=r"$\mathrm{1D}-90\% \chi_\mathrm{st}^{extinction}$") 
ax.plot(zout, rhoY_H2_Z / rho_Z, color="r", linewidth=3.0, linestyle = "-",
        label=r"$\mathrm{3D}$")
ax.plot(zout, states_0D_min.Y[:,H2_ID], color="k", linestyle="--", label="Pure mix")
ax.plot(zout, states_0D_max.Y[:,H2_ID], color="k", linewidth=3.0, linestyle = "-", 
        label="Equilibirum")
ax.legend()
ax.legend(fontsize=20)
ax.set_xlim([0, 0.1])
ax.set_ylim([0, 0.1])
ax.set_ylabel(r"$\langle Y_\mathrm{H2} |Z \rangle$", fontsize=20)
ax.set_xlabel(r"$\langle Z \rangle$", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)

#%%

axis_sum = (1, 2)
wt_xZ_sum = np.sum(wt_xyzZ_sum, axis=axis_sum)
T_xZ_sum = np.sum(T_xyzZ_wtsum, axis=axis_sum)
mf_xZ_sum = np.sum(mf_xyzZ_wtsum, axis=axis_sum)
hrr_xZ_sum = np.sum(hrr_xyzZ_wtsum, axis=axis_sum)

T_Z = T_Z_sum / wt_Z_sum 
mf_Z = mf_Z_sum / wt_Z_sum 
hrr_Z = hrr_Z_sum / wt_Z_sum 

T_xZ = T_xZ_sum / wt_xZ_sum 
mf_xZ = mf_xZ_sum / wt_xZ_sum 
hrr_xZ = hrr_xZ_sum / wt_xZ_sum 

fig, ax = plt.subplots()
ax.plot(zout, T_xZ[0,:], label=r"$x\in [0, 0.25] L_x$")
ax.plot(zout, T_xZ[1,:], label=r"$x\in [0.25, 0.5] L_x$")
ax.plot(zout, T_xZ[2,:], label=r"$x\in [0.5, 0.75] L_x$")
ax.plot(zout, T_xZ[3,:], label=r"$x\in [0.75, 1.0] L_x$")

ax.set_ylim([300, 2500])
ax.set_ylabel(r"$\langle T|Z \rangle$", fontsize=20)

fig, ax = plt.subplots()
ax.plot(zout, hrr_xZ[0,:], label=r"$x\in [0, 0.25] L_x$")
ax.plot(zout, hrr_xZ[1,:], label=r"$x\in [0.25, 0.5] L_x$")
ax.plot(zout, hrr_xZ[2,:], label=r"$x\in [0.5, 0.75] L_x$")
ax.plot(zout, hrr_xZ[3,:], label=r"$x\in [0.75, 1.0] L_x$")
ax.legend(fontsize=20)
#ax.set_ylim([300, 2500])
ax.set_ylabel(r"$\langle HRR|Z \rangle$", fontsize=20)

#%%
fig, ax = plt.subplots()
ax.plot(zout, wt_xZ_sum[0,:]/np.sum(wt_xZ_sum[0,:]), label=r"$x\in [0, 0.25] L_x$")
ax.plot(zout, wt_xZ_sum[1,:]/np.sum(wt_xZ_sum[1,:]), label=r"$x\in [0, 0.25] L_x$")
ax.plot(zout, wt_xZ_sum[2,:]/np.sum(wt_xZ_sum[2,:]), label=r"$x\in [0, 0.25] L_x$")
ax.plot(zout, wt_xZ_sum[3,:]/np.sum(wt_xZ_sum[3,:]), label=r"$x\in [0, 0.25] L_x$")
ax.set_ylabel(r"$PDF(Z)$")
ax.set_ylim([0, 0.1]) 


#%%

#ax.set_

#T = T_wtsum / wt_sum
#%%
dtt = np.array([[9315, 1.118897468e-08],
                [9324, 9.45E-09],
                [9400, 8.409409194e-09],
                [9419, 9.083199516e-09],
                [9450, 5.744984252e-09],
                [9500, 4.819243912e-09],
                [9528, 4.520208237e-09],
                [9567, 6.367467912e-09],
                [9593, 5.938937554e-09],
                [9600, 5.439851261e-09],
                [9615, 6.07375192e-09],
                [9652, 5.259508211e-09],
                [9620, 5.711556162e-09],
                [9700, 4.986878851e-09],
                [9800, 3.882772552e-09],
                [9900, 3.852681304e-09 ],
                [10000, 4.26972249e-09],
                [10093, 5.284692396e-09],
                [10200, 6.722227195e-09],
                [10300, 2.124851619e-08],
                [10400, 2.108682278e-08],
                [10500, 2.030826142e-08],
                [10560, 1.845947361e-08],
                [10570, 2.030826142e-08],
                [10580, 8.214915191e-09],
                [10600, 5.50858402e-09],
                [10700, 5.233702496e-09],
                [10800, 4.447608668e-09],
                [10900, 1.206267925e-09],
                [11040, 8.748198511e-10],
                ])
fig, ax = plt.subplots()
ax.plot(dtt[:,0], dtt[:,1], "bo-")
#%%
#
A = 0.4
B = 0.28
J = 6.66
Djet = 4.5E-4
xs = np.linspace(0, 0.01, 1000)
ys = J * Djet * A * np.power(xs / (J*Djet), B)
#%%
# Declaration
wt_sum = np.zeros((nx, ny, nz))
rho_wtsum = np.zeros((nx, ny, nz))
T_wtsum = np.zeros((nx, ny, nz))
HRR_wtsum = np.zeros((nx, ny, nz))
ts11_wtsum = np.zeros((nx, ny, nz))
ts22_wtsum = np.zeros((nx, ny, nz))
ts33_wtsum = np.zeros((nx, ny, nz))
ts12_wtsum = np.zeros((nx, ny, nz))
ts13_wtsum = np.zeros((nx, ny, nz))
ts23_wtsum = np.zeros((nx, ny, nz))
mu_wtsum = np.zeros((nx, ny, nz))
ts_wtsum = np.zeros((nx, ny, nz))

for ifn, fn in enumerate(fns[0:]):
  f = h5py.File(fn, 'r+')
  wt_sum = wt_sum + np.array(f["DATA"]["volume_mean"])
  rho_wtsum = rho_wtsum + np.array(f["DATA"]["rho_mean"]) 
  T_wtsum = T_wtsum + np.array(f["DATA"]["temp_mean"]) 
  HRR_wtsum = HRR_wtsum + np.array(f["DATA"]["HeatRelease_mean"]) 
  ts11_wtsum = ts11_wtsum + np.array(f["DATA"]["ts11_mean"]) 
  ts22_wtsum = ts22_wtsum + np.array(f["DATA"]["ts22_mean"]) 
  ts33_wtsum = ts33_wtsum + np.array(f["DATA"]["ts33_mean"]) 
  ts12_wtsum = ts12_wtsum + np.array(f["DATA"]["ts12_mean"]) 
  ts13_wtsum = ts13_wtsum + np.array(f["DATA"]["ts13_mean"]) 
  ts23_wtsum = ts23_wtsum + np.array(f["DATA"]["ts23_mean"]) 
  mu_wtsum = mu_wtsum + np.array(f["DATA"]["mu_mean"]) 
  ts_wtsum = ts_wtsum + np.array(f["DATA"]["ts_sum_mean"]) 

rho = rho_wtsum / wt_sum
T = T_wtsum / wt_sum
HRR = HRR_wtsum / wt_sum
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

# T - y
fig, ax = plt.subplots()
vmin = 500; vmax = 2800
im = ax.imshow(T[:,32,:].transpose(), cmap="jet", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
ax.plot(xs, ys)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$ T \; [\mathrm{K}]$", fontsize=20)
plt.savefig("T.png", dpi=300, bbox_inches="tight")
#%%
# HRR - y
fig, ax = plt.subplots()
vmin = 0; vmax = 1E11
phi = HRR[:,27,:]
im = ax.imshow(phi.transpose(), cmap="hot", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("$HRR \; [\mathrm{K}]$", fontsize=20)
plt.savefig("T.png", dpi=300, bbox_inches="tight")

#%%
fig, ax = plt.subplots()
vmin = 0; vmax = 5
im = ax.imshow(rho[:,8,:].transpose(), cmap="jet", extent=[xmin,xmax,zmin,zmax],
               origin="lower",
               vmin = vmin, vmax = vmax,)
#               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
plt.savefig("eta.png", dpi=300, bbox_inches="tight")#
#%%