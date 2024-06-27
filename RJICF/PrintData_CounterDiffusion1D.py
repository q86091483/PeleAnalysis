#%% Loading modules
import sys
import os
import os.path as path
os.environ['MPLCONFIGDIR'] = "./tmp"

sys.path.append(".")
from mixture_fraction import mf
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import cantera as ct; ct.suppress_thermo_warnings()
import graphviz

#%% Read input
# What tasks to do
#do_things = ["do_plot", "do_laminar1D"]
do_things = ["do_laminar1D", "do_plot"]
# Which mixture to calculate
Zs = np.linspace(0, 0.3, 61)
#Zs = np.array([0.0252])
# Which mixture to plot
plot_Zs = np.array([2.483E-01])
# Which mechanism to use
mech = "nuig_H2_4atm"
# Folders
# The relative and absolute path to PeleAnalysis
path_PeleAnalysis = "../.."
abspath_PeleAnalysis = path.abspath(path.join(__file__, path_PeleAnalysis))
# Output folder
output_folder = "Data/CounterDiffusion1D/"
if not os.path.exists(output_folder):
  os.mkdir(output_folder)
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
def get_flame_name(strainrate):
  str1 = "a=" + "%.2E"% strainrate
  str1 = path.join(output_folder, str1)
  str1 = str1 + ".csv"
  return str1

#%% Initialize gas objects and parameters
param = get_param_base()
gas_mix = ct.Solution(param["mech"])
gas_j = ct.Solution(param["mech"])
gas_c = ct.Solution(param["mech"])
species_names = gas_j.species_names
Ns     = len(gas_j.species_names)
i_N2   = gas_j.species_index("N2") 
i_NO   = gas_j.species_index("NO") 
i_NO2  = gas_j.species_index("NO2") 
i_N2O  = gas_j.species_index("N2O") 
i_NNH  = gas_j.species_index("NNH") 
X_j    = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"] 
X_c    = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79
gas_j.TPX = param["T_j"], param["P"], param["X_j"]
gas_c.TPX = param["T_c"], param["P"], param["X_c"]
id_H2 = gas_j.species_index("H2")
id_O2 = gas_j.species_index("O2")

#%% Calculate 1D unstrained flames
param_base = get_param_base()
if "do_laminar1D" in do_things:
  gas = ct.Solution(param["mech"])
  width = 1.0E-3
  f = ct.CounterflowDiffusionFlame(gas, width=width)
  f.P = param_base["P"]  # 1 bar
  f.fuel_inlet.mdot = 20.9  # kg/m^2/s
  f.fuel_inlet.X = 'H2:1'
  f.fuel_inlet.T = 300  # K
  f.oxidizer_inlet.mdot = 125.4 # kg/m^2/s
  f.oxidizer_inlet.X = 'O2:0.21,N2:0.79'
  f.oxidizer_inlet.T = 750  # K

  # Set refinement parameters, if used
  f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, )
  f.solve(loglevel = 0, auto = True)
  f.transport_model = "multicomponent" 
  f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, )
  f.solve(loglevel = 0, auto = True)
  #f.soret_enabled = True
  #f.solve(loglevel = 0, auto = True)


  # Define a limit for the maximum temperature below which the flame is
  # considered as extinguished and the computation is aborted
  # This increases the speed of refinement, if enabled
  temperature_limit_extinction = 900 
  strainrate = np.abs(f.velocity[0]-f.velocity[-1]) / width
  print("strainrate: ", strainrate)
  print("flame name:", get_flame_name(strainrate))
  f.save(get_flame_name(strainrate), basis="mass", overwrite=True)

#%% Plot selected flames
if "do_plot" in do_things:

  #strainrate = 2.09E4;
  fn = get_flame_name(strainrate)
  f = pd.read_csv(fn, index_col=False)
  df = pd.read_csv(fn)
  fstate = ct.SolutionArray(gas_mix)
  fstate.from_pandas(df)
  hrr = -np.sum(fstate.net_production_rates*fstate.partial_molar_enthalpies, axis=1)
  mixfrac = [] 

  gas.X = "O2:0.21, N2:0.79"
  Yox = gas.Y
  gas.X = "H2:1.0"
  Yfu = gas.Y
  zer = mf(gas, Yox, Yfu)
  for ix, x in enumerate(fstate.grid):
    gas.TPY = fstate.T[ix], fstate.P[ix], fstate.Y[ix,:]
    #mixfrac.append(gas.mixture_fraction(fuel = "H2:1",
    #                                    oxidizer = "O2:0.233, N2:0.768",
    #                                   basis="mass",
    #                                   element="Bilger"))
    mixfrac.append(zer.spec2mf(gas.Y))
  fig, ax = plt.subplots()
  ax.plot(mixfrac, hrr)
  ax.set_ylim([0, 2E11])
  fig, ax = plt.subplots()
  ax.plot(mixfrac, fstate.T)


        
# %%
