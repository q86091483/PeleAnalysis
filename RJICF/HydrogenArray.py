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
do_things = ["Lz_J", "do_laminar1D"]
do_things = []
mech = "nuig_H2"

#if (False):
#  fig, ax = plt.subplots()
#  x = [0.813, 0.75, 0.7, 0.619, 0.5935]
#  y = [7.616, 6.46, 5.55, 4.409, 4.05]
#  ax.plot(x, y, color="r", marker="o")
#  ax.set_xlabel(r"$\phi$", fontsize=20)
#  ax.set_ylabel(r"$J$", fontsize=20)
#  ax.tick_params(axis='both', which='major', labelsize=18)
#  ax.tick_params(axis='both', which='minor', labelsize=18)
#  ax.grid()
#%%
# 0. Initialize parameter ---------------------------------------
def get_param_base():
  param_base = {}
  param_base["Re_j"]   = 4913
  param_base["njet"]   = 2
  param_base["P"]      = 10.0 * ct.one_atm 
  param_base["D_j"]    = 5.0E-4
  param_base["J"]      = 6.634717274460312
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
#%%

def get_param_derived(param):
  gas_j = ct.Solution(param["mech"])
  gas_c = ct.Solution(param["mech"])

  id_H2 = gas_j.species_index("H2")
  id_O2 = gas_j.species_index("O2")
  id_N2 = gas_j.species_index("N2")

  gas_j.TPX = param["T_j"], param["P"], param["X_j"]
  gas_c.TPX = param["T_c"], param["P"], param["X_c"]

  rho_j = gas_j.density
  rho_c = gas_c.density

  nu_j = gas_j.viscosity / gas_j.density
  nu_c = gas_c.viscosity / gas_c.density

  U_j = param["Re_j"] * nu_j / param["D_j"]
  A_j = 0.25 * np.pi * param["D_j"] * param["D_j"]
  m_j = param["njet"] * A_j * U_j * rho_j

  Xmix = np.zeros_like(gas_j.X)
  Xmix[id_O2] = 1 / (2 * param["equiv"] + 4.76)
  Xmix[id_H2] = 2 * param["equiv"] * Xmix[id_O2]
  Xmix[id_N2] = 3.76 * Xmix[id_O2]
  gas_mix = ct.Solution(param["mech"])
  gas_mix.X = Xmix
  Ymix = gas_mix.Y
  hmix = gas_j.enthalpy_mass * Ymix[id_H2] + gas_c.enthalpy_mass * (1-Ymix[id_H2])
  gas_mix.HPY = hmix, param["P"], Ymix
  #J = (rho_j * U_j**2) / (rho_c * U_c**2)
  U_c = np.sqrt(1/param["J"]) * np.sqrt(rho_j/rho_c) * U_j
  m_c = m_j * (gas_mix.Y[id_O2] + gas_mix.Y[id_N2])/ gas_mix.Y[id_H2]
  A = m_c / U_c / rho_c
  Ly  = (2*param["intv"]+2) * param["D_j"] 
  Lz = A / Ly
  ufix = (m_j + m_c) / (gas_mix.density * A)
  #print(gas_mix.density, gas_j.density, gas_c.density)
  res = {}
  for k in param.keys():
    res[k] = param[k]
  res["U_j"] = U_j
  res["U_c"] = U_c
  res["Ma_j"] = U_j / gas_j.sound_speed
  res["Ma_c"] = U_c / gas_c.sound_speed
  res["rho_j"] = rho_j
  res["rho_c"] = rho_c
  res["m_j"] = m_j
  res["m_c"] = m_c
  res["ufix"] = ufix
  res["nu_j"] = nu_j
  res["nu_c"] = nu_c
  res["m_j"] = m_j
  res["m_c"] = m_c
  res["Lx"] = 28.4 * param["D_j"]
  res["Ly"] = Ly
  res["Lz"] = Lz
  res["Re_c"] = U_c * Lz / nu_c
  res["h_j"] = gas_j.enthalpy_mass
  res["h_c"] = gas_c.enthalpy_mass
  res["Y_j"] = gas_j.Y
  res["Y_c"] = gas_c.Y

  delta_nu0 = 8.8E-6
  nu0 = nu_c
  u0 = 17.5
  H0 = 6E-3
  alpha = np.sqrt(nu0/(u0/H0))
  alpha = np.sqrt(res["nu_c"]/(res["U_c"]/res["Lz"])) / alpha
  res["delta_nu"] = alpha * delta_nu0
  return res

param_base = get_param_base()
param_derived = get_param_derived(param_base)
#%%
str1 = "Lz_J"
if str1 in do_things:
  intvs = [2, 3, 4, 5]
  Js = [1, 2, 4, 6, 8, 10, 12, 14]
  Lzs = []
  for iintv, intv in enumerate(intvs):
    Lzs.append([])
    for iJ, J in enumerate(Js):
      param_base = get_param_base()
      param_base["J"] = J
      param_base["intv"] = intv
      res = get_param_derived(param_base)
      Lzs[iintv].append(res["Lz"]/res["D_j"])
  fig, ax = plt.subplots()
  for iintv, intv in enumerate(intvs):
    ax.plot(Js, Lzs[iintv], label=r"$L_{intv}=" + "%.1f"%intv + "D_j$")
  ax.set_ylim([0, 25]); ax.set_ylabel(r"$L_z/D_j$", fontsize=20)
  ax.set_xlim([0, 14]); ax.set_xlabel(r"$J$", fontsize=20)
  ax.legend(fontsize=18); ax.grid()
  ax.tick_params(axis='both', which='major', labelsize=16)
  ax.tick_params(axis='both', which='minor', labelsize=16)
  plt.savefig("Lzs.png", dpi=300, bbox_inches="tight")

def get_ngrid(lev, nx0, ny0, nz0, r):
  n = nx0 * ny0 * nz0 * r * np.power(8, lev)
  return n
#%%
def get_flame_name(Patm, Tj, Tc, z):
  str1 = "laminar1D/" + mech + "/P=" + "%i" % Patm
  str1 = str1 + "_Tj=" + "%i" % Tj
  str1 = str1 + "_Tc=" + "%i" % Tc
  str1 = str1 + "_Z=" + "%.4f" % z
  str1 = str1 + ".csv"
  return str1
#%%
Patms = [3, 4, 5, 6, 8, 10]
Tjs = [300, 350]
Tcs = [750]
Zs = [0.0252]
do_laminar1D = False
if "do_laminar1D" in do_things:
  do_laminar1D = True
if do_laminar1D:
  for iP, Patm in enumerate(Patms):
    for m, Tj in enumerate(Tjs):
      for n, Tc in enumerate(Tcs): 
        for k, Z in enumerate(Zs):
          param_base = get_param_base()
          param_base["P"] = Patm * ct.one_atm
          param_base["T_j"] = Tj
          param_base["T_c"] = Tc
          param = get_param_derived(param_base)
          #for k in param.keys():
          #  print(k, ":", param[k])

          hmix = param["h_j"] * Z + param["h_c"] * (1-Z)
          Ymix = param["Y_j"] * Z + param["Y_c"] * (1-Z)
          gas_zst = ct.Solution(param_base["mech"])
          gas_zst.HPY = hmix, Patm*ct.one_atm, Ymix
          f = ct.FreeFlame(gas_zst)
          loglevel = 0
          f.set_max_grid_points(domain=f.domains[1], npmax=10000)
          f.transport_model = 'mixture-averaged'
          f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
          f.solve(loglevel, auto=True)
          f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)
          f.solve(loglevel)  # don't use 'auto' on subsequent solves

          gradf = np.zeros_like(f.T)
          gradf[1:-1] = (f.T[2:]-f.T[0:-2]) / (f.flame.grid[2:]-f.flame.grid[0:-2])
          max_gradf = np.amax(gradf)
          lf = (f.T[-1]-f.T[0]) / max_gradf
          #print("Flame thickness: ", Z, lf)
          #print("Flame speed: ", Z, f.velocity[0])
          #print("nu:", param_derived["nu_j"])

          fn = get_flame_name(Patm, Tj, Tc, Z)
          f.save(fn, basis="mole", overwrite=True)
          #print(fn)
#%% Parameters for computation cost
cost_per_cell_step = 3.0E-7
rs = np.array([np.mean([1.0]), 
              np.mean([1.0]), 
              np.mean([0.8]), 
              np.mean([0.4]),
              np.mean([0.2])  ])
def get_ncell_lev(lev, nx0, ny0, nz0, r):
  n = nx0 * ny0 * nz0 * r * np.power(8, lev)
  return n
def get_total_cell(finest_level, nx0, ny0, nz0):
  num_total_cell = 0
  for lev in range(0, finest_level+1):
    r = rs[lev]
    num_total_cell += get_ncell_lev(lev, nx0, ny0, nz0, r)
  return num_total_cell
#%%
reds=["lightcoral", "red", "brown"]
blues=["lightskyblue", "blue", "midnightblue"]
do_plot_all = True
Tj = 300
Tc = 750
Z = 0.0252
Patms = [3, 4, 5, 6, 8, 10]
Reys = [4000, 4500, 4750]
D_js = [4.5E-4, 5E-4]
ngrid_per_Dj = 5*9/10 # grids for 1 jet diameter
eta=9E-6

lfs = []
sls = []
nu_js = []
U_js = []
U_cs = []
Ma_js = []
Ma_cs = []
ngs3 = []
ngs4 = []
Lzs = []
dts3 = []
dts4 = []
Nts3 = []
Nts4 = []
ncell3=[]
ncell4=[]
costs3 = []
costs4 = []
neta4 = []
ufix = []
delta_nu = []
nnu4 = []
if do_plot_all:
  for iR, Re in enumerate(Reys):
    lfs.append([])
    sls.append([])
    nu_js.append([])
    U_js.append([])
    U_cs.append([])
    Ma_js.append([])
    Ma_cs.append([])
    ngs3.append([])
    ngs4.append([])
    Lzs.append([])
    dts3.append([])
    dts4.append([])
    Nts3.append([])
    Nts4.append([])
    ncell3.append([])
    ncell4.append([])
    costs3.append([])
    costs4.append([])
    neta4.append([])
    ufix.append([])
    delta_nu.append([])
    nnu4.append([])
    for iD, D_j in enumerate(D_js):
      lfs[iR].append([])
      sls[iR].append([])
      nu_js[iR].append([])
      U_js[iR].append([])
      U_cs[iR].append([])
      Ma_js[iR].append([])
      Ma_cs[iR].append([])
      ngs3[iR].append([])
      ngs4[iR].append([])
      Lzs[iR].append([])
      dts3[iR].append([])
      dts4[iR].append([])
      Nts3[iR].append([])
      Nts4[iR].append([])
      ncell3[iR].append([])
      ncell4[iR].append([])
      costs3[iR].append([])
      costs4[iR].append([])
      neta4[iR].append([])
      ufix[iR].append([])
      delta_nu[iR].append([])
      nnu4[iR].append([])
      for iP, Patm in enumerate(Patms):
        # laminar 1D flame
        fn = get_flame_name(Patm, Tj, Tc, Z) 
        f = pd.read_csv(fn)
        T = f["T"].values
        x = f["grid"].values
        u = f["velocity"].values
        gradf = np.zeros_like(T)
        gradf[1:-1] = (T[2:]-T[0:-2]) / (x[2:]-x[0:-2])
        max_gradf = np.amax(gradf)
        #print(fn, max_gradf)
        lf = (T[-1]-T[0]) / max_gradf
        sl = u[0]
        sls[iR][iD].append(sl)
        lfs[iR][iD].append(lf)
        # Physics
        param_base = get_param_base()
        param_base["P"] = Patm * ct.one_atm
        param_base["Re_j"] = Re
        param_base["D_j"] = D_j
        param_derived = get_param_derived(param_base)
        nu_js[iR][iD].append(param_derived["nu_j"])
        U_js[iR][iD].append(param_derived["U_j"])
        U_cs[iR][iD].append(param_derived["U_c"])
        Ma_js[iR][iD].append(param_derived["Ma_j"])
        Ma_cs[iR][iD].append(param_derived["Ma_c"])
        #print(Patm, Re, param_derived["U_j"])
        # Base level mesh size 
        dx0 = param_derived["D_j"] / ngrid_per_Dj
        ngs3[iR][iD].append(lf / (dx0/8))
        ngs4[iR][iD].append(lf / (dx0/16))
        # Domain height
        Lzs[iR][iD].append(param_derived["Lz"]/param_derived["D_j"])
        # Time step size
        cfl = 0.8
        dt3 = cfl * (dx0/8) / param_derived["U_j"]
        dts3[iR][iD].append(dt3)
        dt4 = cfl * (dx0/16) / param_derived["U_j"]
        dts4[iR][iD].append(dt4)
        # Number of time step
        tflow = param_derived["Lx"] / param_derived["U_c"]
        Nt3 = tflow / dt3
        Nt4 = tflow / dt4
        Nts3[iR][iD].append(Nt3)
        Nts4[iR][iD].append(Nt4)
        # Cell number
        ny0 = (param_derived["Ly"]/param_derived["D_j"]) * ngrid_per_Dj
        nx0 = (param_derived["Lx"]/param_derived["Ly"]) * ny0
        nz0 = (param_derived["Lz"]/param_derived["Ly"]) * ny0
        n3 = get_total_cell(3, nx0, ny0, nz0)
        n4 = get_total_cell(4, nx0, ny0, nz0)
        ncell3[iR][iD].append(n3)
        ncell4[iR][iD].append(n4)
        # Cost
        cost3 = cost_per_cell_step * n3 * Nt3
        cost4 = cost_per_cell_step * n4 * Nt4
        costs3[iR][iD].append(cost3)
        costs4[iR][iD].append(cost4)
        # Grid number per Kolmolgorov scale
        eta_Re = eta * np.power(param_derived["Re_j"]/4913 ,-3/4)
        neta4[iR][iD].append(eta_Re / (dx0 / 16))
        # Boundary layer
        ufix[iR][iD].append(param_derived["ufix"])
        #print(param_derived["ufix"])
        delta_nu[iR][iD].append(param_derived["delta_nu"])
        nnu4[iR][iD].append(param_derived["delta_nu"]/(dx0/16))
#%%
field_names = [["lf", "sl", "U_j", "Ma_j", "U_c"], 
               ["ngs", "Lzs", "dts", "Nts", "ncell"],
               ["cost", "neta", "nnu", "ufix", "delta_nu"]]
npy = len(field_names)
npx = len(field_names[0])
fig_unit = 3; figsize = [fig_unit*npx*1.3, fig_unit*npy]
fig, axs = plt.subplots(npy, npx, figsize=figsize)

for ipx in range(0, npx):
  for ipy in range(0, npy):
    var_name = field_names[ipy][ipx]
    ax = axs[ipy, ipx]
    if (var_name == "lf"):
      ax.plot(Patms, np.array(lfs[0][0])/1E-6, linestyle="-", color="k", marker='.')
      ax.set_title(r"$\delta_f \; \mathrm{[\mu m]}$", fontsize=18)
      ax.set_ylim([0.0, 1.5E-4])
      ax.set_yticks([0, 50, 100, 150])
      #ax.set_yticklabels(["0.0", "%.1e"%(5E-5), "%.1e"%(1E-4), "%.1e"%(1.5E-4)])
    elif(var_name == "sl"):
      ax.plot(Patms, sls[0][0], color="k", marker='.')
      ax.set_title(r"$S_l \; \mathrm{[m/s]}$", fontsize=18)
    elif(var_name == "nu_j"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, nu_js[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, nu_js[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$\nu_j \; \mathrm{[m^2/s]}$", fontsize=18)
    elif(var_name == "U_j"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, U_js[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, U_js[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$U_j \; \mathrm{[m/s]}$", fontsize=18)
      ax.set_ylim([0, 400])
    elif(var_name == "U_c"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, U_cs[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, U_cs[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$U_c \; \mathrm{[m/s]}$", fontsize=18)
      ax.set_ylim([0, 60])
    elif(var_name == "Ma_j"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, Ma_js[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, Ma_js[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$Ma_j$", fontsize=18)
      ax.set_ylim([0, 0.3])
    elif(var_name == "ngs"):
      #ax.plot(Patms, ngs3[iR][0], color="r", linestyle="-", marker='.', label=r"$D_j="+"%.1f"%(D_js[0]/1E-3)+"-lev3$")
      #ax.plot(Patms, ngs3[iR][1], color="b", linestyle="-", marker='.', label=r"$D_j="+"%.1f"%(D_js[1]/1E-3)+"-lev3$")
      ax.plot(Patms, ngs4[iR][0], color="r", linestyle="-", marker='.',  label=r"$D_j="+"%.1f"%(D_js[0]/1E-3)+"-lev4$")
      ax.plot(Patms, ngs4[iR][1], color="b", linestyle="-", marker='.',  label=r"$D_j="+"%.1f"%(D_js[1]/1E-3)+"-lev4$")
      ax.set_title(r"$Grid \; \# per \;\delta_f$", fontsize=18)
      ax.set_ylim([0, 30])
      ax.legend(fontsize=8)
    elif(var_name == "Lzs"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, Lzs[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, Lzs[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$L_z/D_j$", fontsize=18)
      ax.set_ylim([11, 13])
    elif(var_name == "dts"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, dts4[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, dts4[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$dt$", fontsize=18)
      #ax.set_ylim([11, 13])
    elif(var_name == "Nts"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, Nts4[iR][0], color=reds[iR], marker='.')
        ax.plot(Patms, Nts4[iR][1], color=blues[iR], marker='.')
      ax.set_title(r"$N_t$", fontsize=18)
      #ax.set_ylim([11, 13])
    elif(var_name == "ncell"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(ncell4[iR][0])/1E6, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(ncell4[iR][1])/1E6, color=blues[iR], marker='.')
      ax.set_title(r"$Total \; cell \; number \; \mathrm{[M]}$", fontsize=18)
    elif(var_name == "cost"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(costs4[iR][0])/1E6, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(costs4[iR][1])/1E6, color=blues[iR], marker='.')
      ax.set_title(r"$Cost \; in \; CPUh \; \mathrm{[M]}$", fontsize=18)
      ax.set_ylim([0, 2.5])
    elif(var_name == "neta"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(neta4[iR][0])/1.0, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(neta4[iR][1])/1.0, color=blues[iR], marker='.')
      ax.set_ylim([0, 2.5], )
      ax.set_title(r"$Grid \; \# \; per \; \eta_k$", fontsize=18)
    elif(var_name == "ufix"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(ufix[iR][0])/1.0, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(ufix[iR][1])/1.0, color=blues[iR], marker='.')
      ax.set_title(r"$U_{c,fix}$", fontsize=18)
      ax.set_ylim([0, 60])
    elif(var_name == "nnu"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(nnu4[iR][0])/1.0, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(nnu4[iR][1])/1.0, color=blues[iR], marker='.')
      ax.set_title(r"$Grid \; \# \; per \; \delta_\nu$", fontsize=18)
      #ax.set_ylim([0, 2.5])
    elif(var_name == "delta_nu"):
      for iR in range(0, len(Reys)):
        ax.plot(Patms, np.array(delta_nu[iR][0])/1.0E-6, color=reds[iR], marker='.')
        ax.plot(Patms, np.array(delta_nu[iR][1])/1.0E-6, color=blues[iR], marker='.')
      ax.set_title(r"$\delta_\nu \; \mathrm{[\mu m]}$", fontsize=18)
      #ax.set_ylim([0, 60])


    ax.grid(linestyle="--", color="gray")
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    ax.set_xlabel(r"$P \; \mathrm{[atm]}$", fontsize=18)
    ax.set_xticks([4,6,8,10])
    ax.set_xlim([3,10])
    if ipy < npy-1:
      ax.set_xlabel("")
      ax.set_xticklabels("")
plt.savefig("params.png", dpi=300, bbox_inches="tight")

#%%


      

