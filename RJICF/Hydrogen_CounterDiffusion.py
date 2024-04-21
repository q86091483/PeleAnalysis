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
import graphviz

#%%
# Input
# Tasks to perform
do_things = ["Lz_J", "do_laminar1D"]
do_things = []
# Mechanism
Zs = np.array([0.0252])
mech = "nuig_H2_4atm"
# Output folder
output_folder = "CounterDiffusionPremixed1D/"
if not os.path.exists(output_folder):
  os.mkdir(output_folder)
#%%
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

def get_flame_name(z, a):
  str1 = "Z=" + "%.3E" % z
  str1 = str1 + "_a=" + "%.3E" % a
  str1 = str1 + ".csv"
  str1 = os.path.join(output_folder, str1)
  return str1

#%%
param = get_param_base()
gas_mix = ct.Solution(param["mech"])
gas_j = ct.Solution(param["mech"])
gas_c = ct.Solution(param["mech"])

Ns     = len(gas_j.species_names)
id_H2 = gas_mix.species_index("H2")
id_O2 = gas_mix.species_index("O2")
i_N2   = gas_j.species_index("N2") 
i_NO   = gas_j.species_index("NO") 
i_NO2  = gas_j.species_index("NO2") 
i_N2O  = gas_j.species_index("N2O") 
i_NNH  = gas_j.species_index("NNH")

X_j    = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"] 
X_c    = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79

gas_j.TPX = param["T_j"], param["P"], param["X_j"]
gas_c.TPX = param["T_c"], param["P"], param["X_c"]

for iz, z in enumerate(Zs):
  Ys = z*gas_j.Y + (1-z)*gas_c.Y
  Hs = z*gas_j.enthalpy_mass + (1-z)*gas_c.enthalpy_mass
  Ps = gas_j.P
  gas_mix.HPY = Hs, Ps, Ys
  T_mix = gas_mix.T


  # Create the flame simulation object
  width = 1E-2
  U_in = 3.
  mdot = gas_mix.density * U_in
  sim = ct.CounterflowPremixedFlame(gas=gas_mix, width=width)
  sim.reactants.mdot = mdot
  sim.products.mdot = mdot
  sim.set_refine_criteria(ratio=3, slope=0.1, curve=0.2, prune=0.02)
  sim.set_initial_guess()
  sim.show()
  
  loglevel = 1
  sim.solve(loglevel, auto=True)

#%%