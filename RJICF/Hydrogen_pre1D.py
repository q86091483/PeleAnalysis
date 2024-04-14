#%%
import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = "./tmp"
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import cantera as ct
import pandas as pd

#%%
# Input
# Tasks to perform
do_things = ["Lz_J", "do_laminar1D"]
do_things = []
# Mechanism
Zs = np.array([0.0252])
mech = "nuig_H2_4atm"

def get_param_base():
  param_base = {}
  param_base["Re_j"]   = 3200
  param_base["njet"]   = 2
  param_base["P"]      = 4.0 * ct.one_atm 
  param_base["D_j"]    = 4.5E-4
  param_base["J"]      = 6.634
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

def get_flame_name(Patm, z):
  str1 = "laminar1D/" + mech + "/P=" + "%i" % Patm
  str1 = str1 + "_Z=" + "%.4f" % z
  str1 = str1 + ".csv"
  return str1

#%%
param = get_param_base()

gas_mix = ct.Solution(param["mech"])
gas_j = ct.Solution(param["mech"])
gas_c = ct.Solution(param["mech"])

Ns     = len(gas_j.species_names)
X_j    = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"] 
X_c    = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79

gas_j.TPX = param["T_j"], param["P"], param["X_j"]
gas_c.TPX = param["T_c"], param["P"], param["X_c"]

id_H2 = gas_j.species_index("H2")
id_O2 = gas_j.species_index("O2")

fig, axz = plt.subplots(figsize=[6, 4.8])
for iz, z in enumerate(Zs):
  Ys = z*gas_j.Y + (1-z)*gas_c.Y
  Hs = z*gas_j.enthalpy_mass + (1-z)*gas_c.enthalpy_mass
  Ps = gas_j.P
  gas_mix.HPY = Hs, param["P"], Ys
  T_mix = gas_mix.T



  fn = get_flame_name(param["P"]/ct.one_atm, z)
  if "do_laminar1D" in do_things:
    f = ct.FreeFlame(gas_mix)
    loglevel = 0
    f.set_max_grid_points(domain=f.domains[1], npmax=10000)
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel)  # don't use 'auto' on subsequent solves
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)
    f.save(fn, basis="mass", overwrite=True)
  else:
    f = pd.read_csv(fn, index_col=False)

  x_f = f["grid"].values * 1000
  T_f = f["T"].values
  Y_f = []
  for isp, spn in enumerate(gas_mix.species_names):
    Y_f.append(f["Y_"+spn].values)
  grad_T = np.gradient(T_f, x_f)
  i_Tmax = np.argmax(grad_T) 

  states = ct.SolutionArray(gas_mix)
  T_s = T_f[i_Tmax]
  Y_s = {}
  for spn in gas_mix.species_names:
    Y_s[spn] = Y_f[isp][i_Tmax]
  P_s = gas_mix.P
  states.append(T=T_s,
                P=P_s,
                Y=Y_s)
  rates = states.net_rates_of_progress
  

  vmin=700; vmax=2500
  axz.plot(x_f, T_f, label="T", linestyle="--", color="r")
  axz.plot([x_f[i_Tmax], x_f[i_Tmax]], [vmin, vmax], "k--")
  axz.set_xlim([x_f[-1]*0.298, x_f[-1]*0.305])
  axz.set_ylim([vmin, vmax])
  axz.set_xlabel(r"$x \; [mm]$", fontsize=20)
  axz.set_ylabel(r"$T \; [K]$", fontsize=20)
  axz.legend(fontsize=14)

  axz2 = axz.twinx()
  spn = "NNH"; isp = gas_mix.species_index(spn) 
  axz2.plot(x_f, Y_f[isp]*(10**4), label =r"$Y_\mathrm{"+spn+"} * 10^{4}$")
  spn = "NO"; isp = gas_mix.species_index(spn) 
  axz2.plot(x_f, Y_f[isp]*(10**2), label =r"$Y_\mathrm{"+spn+"}*10^{2}$")
  spn = "N2O"; isp = gas_mix.species_index(spn) 
  axz2.plot(x_f, Y_f[isp]*(2*10**3), label =r"$Y_\mathrm{"+spn+"} * 2 * 10^{3}$")
  spn = "OH"; isp = gas_mix.species_index(spn) 
  axz2.plot(x_f, Y_f[isp], label =r"$Y_\mathrm{"+spn+"}$")


  vmin = 600; vmax = 2500
  axz2.set_ylim([0, 2E-2])
  axz2.legend(fontsize=14)



#%%
