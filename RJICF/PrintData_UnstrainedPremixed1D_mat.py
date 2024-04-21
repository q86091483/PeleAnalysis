#%% Loading modules
import sys
import os
import os.path as path
os.environ['MPLCONFIGDIR'] = "./tmp"

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import re
import glob

import cantera as ct; ct.suppress_thermo_warnings()
import graphviz

#%% Read input
# What tasks to do
#do_things = ["do_matrix"]
do_things = ["do_matrix"]
# Which flame.csv to process
csv_pattern = "Z=*.csv"
# Which mechanism to use
mech = "nuig_H2_4atm"
# Folders
# The relative and absolute path to PeleAnalysis
path_PeleAnalysis = "../.."
abspath_PeleAnalysis = path.abspath(path.join(__file__, path_PeleAnalysis))
# Data folder - where flame data is stored
data_folder = "Data/UnstrainedPremixed1D/"
# Output folder
fig_folder = "Figure/UnstrainedPremixed1D/"
if not os.path.exists(fig_folder):
  os.mkdir(fig_folder)
# Files to read
fns_unsorted = glob.glob(path.join(data_folder, csv_pattern))
def get_Z_key(s):
  ss = re.split("/", s)[-1]
  ss = re.split("=", ss)[-1]
  ss = re.split(".csv", ss)[0]
  return(float(ss))
fns_sorted = sorted(fns_unsorted, key=get_Z_key)
# Record mixture fraction array
Zs = []
for ifn, fn in enumerate(fns_sorted):
  Zs.append(get_Z_key(fn))
# Pararameters of the RJICF case
def get_param_base():
  param_base = {}
  param_base["Re_j"]   = 4000
  param_base["njet"]   = 2
  param_base["P"]      = 4.0 * ct.one_atm 
  param_base["D_j"]    = 4.5E-4
  param_base["J"]      = 6.6
  param_base["equiv"]  = 0.6
  param_base["intv"]   = 3.0
  param_base["T_j"]    = 300.
  param_base["T_c"]    = 750.
  X_j = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"] 
  X_c = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79
  param_base["X_j"]    = X_j
  param_base["X_c"]    = X_c
  param_base["mech"]   = mech + "/chem.yaml"
  return param_base
# Name convention
def get_flame_csv(z):
  str1 = "Z=" + "%.3E"% z
  str1 = path.join(data_folder, str1)
  str1 = str1 + ".csv"
  return str1

#%% Initialize gas objects and parameters
param = get_param_base()
gas_mix = ct.Solution(param["mech"])
gas_j = ct.Solution(param["mech"])
gas_c = ct.Solution(param["mech"])
species_names = gas_j.species_names
Ns     = len(gas_j.species_names)
i_H2O   = gas_j.species_index("H2O") 
i_H2 = gas_j.species_index("H2")
i_O2 = gas_j.species_index("O2")
i_N2   = gas_j.species_index("N2") 
i_NO   = gas_j.species_index("NO") 
i_NO2  = gas_j.species_index("NO2") 
i_N2O  = gas_j.species_index("N2O") 
i_NNH  = gas_j.species_index("NNH") 
X_j    = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"] 
X_c    = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79
gas_j.TPX = param["T_j"], param["P"], param["X_j"]
gas_c.TPX = param["T_c"], param["P"], param["X_c"]

#%% Progress variable
coeff_pv = {}
for isp, spn in enumerate(species_names):
  if ("N" in spn or spn=="H2" or spn=="O2" or spn=="N2" or spn=="OHV"):
    coeff_pv[spn] = 0.0
  else:
    coeff_pv[spn] = 1.0
  #if (spn=="H2O"):
  #  coeff_pv[spn] = 1.0
  #elif (spn=="H2"):
  #  coeff_pv[spn] = -1.0
  #elif (spn=="O2"):
  #  coeff_pv[spn] = -1.0
  #else:
  #  coeff_pv[spn] = 0.0

  print(spn)
 
def get_pveq(zs):
  gas_f = ct.Solution(param["mech"])
  gas_o = ct.Solution(param["mech"])
  gas_m = ct.Solution(param["mech"])
  states = ct.SolutionArray(gas_m)
  print(zs)
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
  states.equilibrate("HP")

  pveq = np.zeros_like(states.T)
  for isp, spn in enumerate(gas_f.species_names):
    pveq = pveq + coeff_pv[spn] * states.Y[:,isp]
  return pveq

#%% Plot selected flames

pveqs = get_pveq(Zs)
#%%
fig, ax = plt.subplots()
fig, axc = plt.subplots()

if "do_matrix" in do_things:
  for ifn, fn in enumerate(fns_sorted):
    z = Zs[ifn]
    f = pd.read_csv(fn, index_col=False)
    df = pd.read_csv(fn)
    fstate = ct.SolutionArray(gas_mix)
    fstate.from_pandas(df)
    hrr = -np.sum(fstate.net_production_rates*fstate.partial_molar_enthalpies, axis=1)
    x_f = fstate.grid * 1000
    T_f = fstate.T
    Y_f = []
    for isp, spn in enumerate(gas_mix.species_names):
      Y_f.append(fstate.Y[:,isp])

    pv = np.zeros_like(T_f)
    for isp, spn in enumerate(species_names):
      pv = pv + coeff_pv[spn] * fstate.Y[:,isp]
    pv_min = np.amin(pv)
    pv_max = np.amax(pv)
    pv = (pv-pv_min) / (pv_max - pv_min)
    #for isp, spn in enumerate(gas_mix.species_names):
    #  pv = pv + coeff_pv[spn] * Y_f[isp] 
    #pv = pv / pveqs[ifn]

    grad_T = np.gradient(T_f, x_f)
    i_Tmax = np.argmax(grad_T) 
    delta_T = (T_f[-1]-T_f[0])/grad_T[i_Tmax]
 
    # Plot - xyz space
    vmin=700; vmax=2800
    ax.plot(x_f, T_f, label="T", linestyle="--", color="r")
    #ax.plot([x_f[i_Tmax], x_f[i_Tmax]], [vmin, vmax], "k--")
    ax.set_xlim(x_f[i_Tmax]-5*delta_T, x_f[i_Tmax]+5*delta_T)
    ax.set_ylim([vmin, vmax])
    ax.set_xlabel(r"$x \; [mm]$", fontsize=20)
    ax.set_ylabel(r"$T \; [K]$", fontsize=20)
    ax.legend(fontsize=14)

    ax2 = ax.twinx()
    ax2.plot(x_f, hrr/np.amax(hrr), label =r"$\mathrm{HRR}$")
    ax2.plot(x_f, pv, label =r"$C$")

    # Plot - C space
    phi = fstate.Y[:,gas_mix.species_index("OH")]
    phi = hrr 
    axc.plot(pv, phi/np.amax(phi))

# %%
