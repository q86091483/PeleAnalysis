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

import cantera as ct; ct.suppress_thermo_warnings()
import graphviz

#%% Read input
# What tasks to do
#do_things = ["do_plot", "do_laminar1D"]
do_things = ["do_laminar1D"]
# Which mixture to calculate
#Zs = np.linspace(0, 0.3, 61)
Zs = np.array([0.0275, 0.025])
# Which mixture to plot
plot_Zs = np.array([2.483E-01])
# Which mechanism to use
mech = "nuig_H2_4atm"
# Folders
# The relative and absolute path to PeleAnalysis
path_PeleAnalysis = "../.."
abspath_PeleAnalysis = path.abspath(path.join(__file__, path_PeleAnalysis))
# Output folder
output_folder = "Data/UnstrainedPremixed1D/"
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
def get_flame_name(z):
  str1 = "Z=" + "%.3E"% z
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
if "do_laminar1D" in do_things:
  for iz, z in enumerate(Zs):
    Ys = z*gas_j.Y + (1-z)*gas_c.Y
    Hs = z*gas_j.enthalpy_mass + (1-z)*gas_c.enthalpy_mass
    Ps = gas_j.P
    gas_mix.HPY = Hs, param["P"], Ys
    T_mix = gas_mix.T

    try:
      f = ct.FreeFlame(gas_mix)
      loglevel = 0
      f.set_max_grid_points(domain=f.domains[1], npmax=10000)
      f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
      f.transport_model = 'mixture-averaged'
      f.solve(loglevel)  
      f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)
      f.solve(loglevel)  
      fn = get_flame_name(z)
      f.save(fn, basis="mass", overwrite=True)
      print("FreeFlame.solve() successed at Z = ", z)
    except:
      print("FreeFlame.solve() failed at Z =", z)

#%% Plot selected flames
if "do_plot" in do_things:
  for iz, z in enumerate(plot_Zs):
    fn = get_flame_name(z)
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
    grad_T = np.gradient(T_f, x_f)
    i_Tmax = np.argmax(grad_T) 

    # Definition of progress variable
    coeff_pv = {}
    idef = 0
    if idef == 0:
      for isp, spn in enumerate(gas_mix.species_names):
        if ("N" in spn or spn=="H2" or spn=="O2" or spn=="N2"):
          coeff_pv[spn] = 0.0
        else:
          coeff_pv[spn] = 1.0
    elif idef == 1:
      for isp, spn in enumerate(gas_mix.species_names):
        if (spn == "H2O"):
          coeff_pv[spn] = 1.0
        elif (spn == "H2"):
          coeff_pv[spn] = -1.0
        elif (spn == "O2"):
          coeff_pv[spn] = -1.0
        else:
          coeff_pv[spn] = 0.0

    # Plot
    fig, ax = plt.subplots()
    vmin=700; vmax=2500
    ax.plot(x_f, T_f, label="T", linestyle="--", color="r")
    #ax.plot([x_f[i_Tmax], x_f[i_Tmax]], [vmin, vmax], "k--")
    #ax.set_xlim([x_f[-1], x_f[-1]])
    ax.set_ylim([vmin, vmax])
    ax.set_xlabel(r"$x \; [mm]$", fontsize=20)
    ax.set_ylabel(r"$T \; [K]$", fontsize=20)
    ax.legend(fontsize=14)

    ax2 = ax.twinx()
    ax2.plot(x_f, hrr/np.amax(hrr), label =r"$\mathrm{HRR}$")

#%%
  #states = ct.SolutionArray(gas_mix)
  #T_s = T_f[i_Tmax]
  #Y_s = {}
  #for isp, spn in enumerate(gas_mix.species_names):
  #  Y_s[spn] = Y_f[isp][i_Tmax]
  #P_s = gas_mix.P
  #states.append(T=T_s,
  #              P=P_s,
  #              Y=Y_s)
  #rates = states.net_rates_of_progress

  # Progress variable
  # Reaction pathway analysis
  #tran_matrix = np.genfromtxt("./nuig_H2_4atm/n-N-sp1-sp2_nuig_H2_32sp_4atm.csv", delimiter=",")
  #graph = np.zeros((Ns, Ns))
  #reactions_of_species = {}
  #for iline in range(0, tran_matrix.shape[0]):
  #  ir    = int(tran_matrix[iline, 0])
  #  isp0  = int(tran_matrix[iline, 1])
  #  isp1  = int(tran_matrix[iline, 2])
  #  nN    = int(tran_matrix[iline, 3])

  #  graph[isp0, isp1] = graph[isp0, isp1] + rates[0][ir] * nN
  #sum = np.sum(graph[i_N2,:])

  #ratio = np.zeros((Ns))
  #for isp, spn in enumerate(gas_mix.species_names):
  #  ratio[isp] = graph[i_N2, isp] - graph[isp, i_N2]
  #  if ratio[isp] < 0:
  #    ratio[isp] = 0.0
  #for isp, spn in enumerate(gas_mix.species_names):
  #  print("N2 -> " + spn, ratio[isp]/np.sum(ratio))

  


  #spn = "NNH"; isp = gas_mix.species_index(spn) 
  #axz2.plot(x_f, Y_f[isp]*(10**4), label =r"$Y_\mathrm{"+spn+"} * 10^{4}$")
  #spn = "NO"; isp = gas_mix.species_index(spn) 
  #axz2.plot(x_f, Y_f[isp]*(10**2), label =r"$Y_\mathrm{"+spn+"}*10^{2}$")
  #spn = "N2O"; isp = gas_mix.species_index(spn) 
  #axz2.plot(x_f, Y_f[isp]*(2*10**3), label =r"$Y_\mathrm{"+spn+"} * 2 * 10^{3}$")
  #spn = "OH"; isp = gas_mix.species_index(spn) 
  #axz2.plot(x_f, Y_f[isp], label =r"$Y_\mathrm{"+spn+"}$")

  #vmin = 600; vmax = 2500
  #axz2.set_ylim([0, 2E-2])
  #axz2.legend(fontsize=14)
#%%
#print(graph)
#fig, ax = plt.subplots()
#f = graphviz.Digraph()
# For every "reactant" specie
#for i in range(0, Ns):
    # For every "product" specie 
#    for j in range(0,Ns):
        
#        spo = gas_mix.species_names[i]
#        origin = f"{spo}"
        
#        spd = gas_mix.species_names[j]
#        destin = f"{spd}"

#        color='black'
#        width = 10*graph[i, j]/graph.max()
#        width = max(1, width)

#        label = f"  {100*graph[i, j]/graph.max():.2f} %\n\n"
#        label = f"   {graph[i, j]:.2e}\n\n"
#        add_edge = (graph[i, j] / graph.max()) > 0.02
#        if add_edge:
#            f.edge(origin, destin, label=label, penwidth=str(width), color=color)
#f.view()
#%%
