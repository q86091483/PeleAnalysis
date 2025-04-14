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

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
labelsize = 24; lw = 1.5

# Input
# Where to output the result
case_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_fluxRPA/"

# Where to read h5 files that contain condition smean
fns = glob.glob(case_folder + "plt_17*.h5")

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
gas1D = ct.Solution(mech)
Nsp = len(gas1D.species_names)
species_names = gas1D.species_names
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

for ifn, fn in enumerate(fns):
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
  for ir in range(0, Nr):
    fdn = "R(" + str(ir) + ")_mean"
    prr_wtsum[:,:,:,:,:,:,:,:,ir] = prr_wtsum[:,:,:,:,:,:,:,:,ir] + f["DATA"][fdn]
  print("Finished: ", fn)
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

#%%
def plot_graph(fg, ax, arr, gas):
  # Plot parameters
  species_names = gas.species_names
  xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0
  labelsize = 14
  bw = 0.1
  bh = 0.08
  #flux_ref = np.sum(np.absolute(arr[N2_ID,:]))
  flux_ref = np.sum(arr[N2_ID,:])

  spns = ["N2", "NO", "N", "NNH", "N2O", "NH", "NH2", "HNO", "HONO",]
  isps = []
  for spn in spns:
    isps.append(gas.species_index(spn))
  loc_sp = {}
  cdn = [0.05, 0.3, 0.6, 0.9]
  loc_sp["N2"]    = np.array((cdn[1], cdn[3]))
  loc_sp["NO"]    = np.array((cdn[1], cdn[2]))
  loc_sp["N"]     = np.array((cdn[0], cdn[2]))
  loc_sp["NNH"]   = np.array((cdn[2], cdn[2]))
  loc_sp["N2O"]   = np.array((cdn[3], cdn[3]))
  loc_sp["NH"]    = np.array((cdn[3], cdn[2]))
  loc_sp["NH2"]   = np.array((cdn[3], cdn[1]))
  loc_sp["HNO"]   = np.array((cdn[1], cdn[1]))
  loc_sp["HONO"]  = np.array((cdn[2], cdn[1]))

  # Plot sub-functions
  def write_label(label_x, label_y, str):
    ax.text(label_x, label_y, str, fontsize = labelsize,
            horizontalalignment='center',
            verticalalignment='center')
    x0 = label_x - bw * 0.5
    x1 = label_x + bw * 0.5
    y0 = label_y - bh * 0.5
    y1 = label_y + bh * 0.5
    ax.plot((x0, x1), (y0, y0), color = "black", linewidth = lw)
    ax.plot((x0, x1), (y1, y1), color = "black", linewidth = lw)
    ax.plot((x0, x0), (y0, y1), color = "black", linewidth = lw)
    ax.plot((x1, x1), (y0, y1), color = "black", linewidth = lw)

  def get_start_end(xy0, xy1, r):
    p0 = np.zeros((2)); p1 = np.zeros((2))
    rx0 = 1.0; rx1 = 1.0; ry0 = 1.0; ry1 = 1.0
    x0 = xy0[0]; y0 = xy0[1]; x1 = xy1[0]; y1 = xy1[1]
    if (r > 0):
      p0[0] = x0; p0[1] = y0
      p1[0] = x1; p1[1] = y1
    else:
      p0[0] = x1; p0[1] = y1
      p1[0] = x0; p1[1] = y0
    if (p0[0] == p1[0]):
      rx0 = 0.0
    elif (p0[0] < p1[0]):
      rx0 = 1.0
    else:
      rx0 = -1.0
    if (p0[1] == p1[1]):
      ry0 = 0.0
    elif (p0[1] < p1[1]):
      ry0 = 1.0
    else:
      ry0 = -1.0
    rx1 = -rx0
    ry1 = -ry0
    p0[0] = p0[0] + rx0 * bw / 2.;
    p1[0] = p1[0] + rx1 * bw / 2.;
    p0[1] = p0[1] + ry0 * bh / 2.;
    p1[1] = p1[1] + ry1 * bh / 2.;
    return (p0, p1)

  def label_flux(p0, p1, isp0, isp1, color):
    ptxt = (p0 + p1) / 2.
    net_rate = arr[isp0, isp1] - arr[isp1, isp0]
    v = (np.abs(net_rate / flux_ref)) * 100.
    txt = "%.2f"%np.abs(v)
    ax.text(ptxt[0], ptxt[1], txt + "%", fontsize = 12, color = color)

  def draw_line(p0, p1, r, color, arrow = False):
    if (arrow):
      #ax.arrow(p0[0], p0[1], p1[0]-p0[0], p1[1]-p0[1], color = color, linewidth = lw, head_width = 0.02)
      ax.annotate("", xy=(p1[0], p1[1]), xytext=(p0[0], p0[1]),
            arrowprops=dict(arrowstyle="->", linewidth = lw, color = color))
    else:
      ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color = color, linewidth = lw)

  def connect(sp0, sp1, r, color, middle_points = []):
    xy0 = loc_sp[sp0]
    xy1 = loc_sp[sp1]
    isp0 = gas.species_index(sp0)
    isp1 = gas.species_index(sp1)
    (p0, p1) = get_start_end(xy0, xy1, r)
    if len(middle_points) == 0:
      draw_line(p0, p1,
              r = r,
              color = color,
              arrow = True)
      label_flux(p0, p1, isp0, isp1, color)
    else:
      if r < 0:
        mps = list(reversed(middle_points))
      else:
        mps = middle_points
      for ip, pm in enumerate(mps):
        p_start = np.zeros((2))
        p_end = pm
        if (ip == 0):
          p_start, p_temp = get_start_end(xy0, pm, r)
          label_flux(xy0, pm, isp0, isp1, color)
        else:
          p_start = mps[ip - 1]
        draw_line(p_start, p_end, r = r, color = color, arrow = False)
      (temp_0, temp_1) = get_start_end(mps[-1], xy1, r)
      draw_line(mps[-1], temp_1, r = r, color = color, arrow = True)

  # Connect
  for spn in spns:
    write_label(loc_sp[spn][0], loc_sp[spn][1], spn)

  s0 = "N2"; s1 = "NO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2"; s1 = "N2O"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2"; s1 = "NNH"; color = "k"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2"; s1 = "N"; color = "k"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NO"; s1 = "N"; color = "k"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NO"; s1 = "NNH"; color = "k"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2O"; s1 = "NNH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2O"; s1 = "NO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N2O"; s1 = "NH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color,
   middle_points = [np.array([cdn[3], cdn[3]])])

  s0 = "NH2"; s1 = "NNH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NH2"; s1 = "NH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NNH"; s1 = "NH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "HNO"; s1 = "HONO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NO"; s1 = "HNO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NO"; s1 = "HONO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "NNH"; s1 = "HONO"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "HONO"; s1 = "NH2"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "HNO"; s1 = "N"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color)

  s0 = "N"; s1 = "NH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color,
   middle_points = [np.array([cdn[0], cdn[0]]), np.array([0.99, cdn[0]]), np.array([0.99, cdn[2]])])

  s0 = "N"; s1 = "NH"; color = "r"
  isp0 = gas.species_index(s0); isp1 = gas.species_index(s1)
  rv = arr[isp0, isp1] - arr[isp1, isp0]
  connect(s0, s1, r = rv, color = color,
  middle_points = [np.array([cdn[3], 0.2]), np.array([cdn[1], 0.2])])

  ax.set_xlim([xmin, xmax])
  ax.set_ylim([ymin, ymax])
  ax.set_xlabel("")
  ax.set_ylabel("")
  ax.set_xticks([])
  ax.set_yticks([])
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)

#%%
# 0 - x
# 1 - y
# 2 - z
# 3 - Z
# 4 - C
# 5 - zone
# 6 - FI
# 7 - tres
axis_all = (0, 1, 2, 3, 4, 5, 6, 7)
axis_0123456 = (0, 1, 2, 3, 4, 5, 6)

# All conditions variables
prr_all = np.zeros((Nr))
for ir in range(0, Nr):
  prr_all[ir] = np.sum(prr_wtsum[:,:,:,:,:,:,:,:,ir],
                       axis = axis_all)
graph_all = get_graph(prr_all)

# All premixed
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

#%%=========== Plot graph ============
fig, ax = plt.subplots(figsize = (6, 6), frameon=False)

gt = get_graph(prr_all)
target_ID = N2_ID
target_name = species_names[target_ID]
for isp in range(0, Nsp):
  print(target_name + " to ".rjust(5)+gas1D.species_names[isp].ljust(5), ": "
         "%5.2E" % (gt[target_ID, isp]) )

plot_graph(fig, ax,
           gt, gas1D)
#%% 1D test
fn1D    = fns1D_sorted[23]
df      = pd.read_csv(fn1D)
fstate  = ct.SolutionArray(gas1D)
fstate.from_pandas(df)
x_1D    = fstate.grid
T_1D    = fstate.T
u_1D    = fstate.velocity
C_1D    = fstate.concentrations
prr_1D  = fstate.net_rates_of_progress
frr_1D  = fstate.forward_rates_of_progress
rrr_1D  = fstate.reverse_rates_of_progress
wdot_1D = fstate.net_production_rates


tID     = NO2_ID
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
#from_others = np.sum(graph_1D[:,tID]) - graph_1D[tID,tID]
#to_others = np.sum(graph_1D[tID,:])
#v3 = (to_others - from_others)/2.
v3 = -(np.sum(graph_1D[tID,:]) - np.sum(graph_1D[:,tID])) / gas1D.n_atoms(gas1D.species_names[tID], "N")

print("1. Outflow flux:", v1, "kmol/m2/s")
print("2. Production of NO:", v2, "kmol/m2/s")
print("3. Production of NO from transfer matrix:", v3, "kmol/m2/s")
fig, ax = plt.subplots()
plot_graph(fig, ax,
           graph_1D, gas1D)

import graphviz

#%%
fg_part = graphviz.Digraph(filename="part_premix_rich",
                           format="png",
                          graph_attr={"nodesep" : "0.0",
                                      "pad" : "0.1",
                                      "fontname" : "times bold"},
)
gt = graph_pre
nodesize = 50
labelsize = 40
flux_ref = np.sum(gt[N2_ID,:])
spg_part = ["N2", "NO", "N", "NNH", "N2O"]
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
    add_edge = (gt[isp0, isp1] / flux_ref) > 0.01
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