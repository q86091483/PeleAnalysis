#%%
import sys
import os
os.environ['MPLCONFIGDIR'] = "./tmp"
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import glob as glob
import NGA_reader as NGA_reader
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
#%%
# Documentation:
# The code convert the NGA output of channel turbulence to inlet of Pele hydrogen JICF.
# There are 2 ways to do the above:
# 1. Frozen input - 1 NGA output (e.g., data_3.000E-02) is periodically fed,
#    like moving the data box with mean cross flow velocity
# 2. Time dependent input - NGA at x_target at different times is fed.

# Inputs:
# Folder that contains the case
case_folder   = "/scratch/b/bsavard/zisen347/cases/"
# Case name
case_name     = "NUIG_Re4000_2J6_4atm"
# Case path to input.3d, turbin of Pele
case_path     = os.path.join(case_folder, case_name)
# The time range t0 to t1 should cover the Pele run
t0            = 5.5E-4    # 2.0E-4
t1            = 1.0E-3    # 6.0E-4
# The NGA output, data_(t_nga_froz), is first periodically fed 
t_nga_froz    = 1.000E-2;
# Time when the channel input of Pele change from frozen flow 
# to time-depedent flow
t_pele_froz2tran = 2.26E-4;
fns = glob.glob(os.path.join(case_path,"data_*"))
# Parameter of JICF
U_c           = 37.5    # [m/s]
Lx            = 1.28E-2 # [m]
# Some input of NGA channel turbulence
dataname = os.path.join(case_path, "data_"+"%.3E"%t_nga_froz)
configname = os.path.join(case_path, "config")
din = NGA_reader.NGA_reader(dataname, configname)
# Find which slice along x direction to output based on t_pele_froz2tran
x_target = din["x"][0] + din["Lx"] - np.fmod(U_c * t_pele_froz2tran, din["Lx"]);
ix_target = np.argmax(din["x"] > x_target)

#%%


# We iterate over time range [t0, t1] + t_nga_froz of NGA channel outputs.
# For each of these time instants, we extract the at x = x_target.
def get_key(fn):
  s = re.split("_", fn)[-1]
  s = float(s)
  return s

fns_sorted = sorted(fns, key=get_key)
ts =[]
for fn in fns_sorted:
  ts.append(get_key(fn))
ts = np.array(ts)

i0 = np.argmax(ts > (t0 + t_nga_froz))
if i0 > 0:
  i0 = i0 - 1
i1 = np.argmax(ts > (t1 + t_nga_froz))

lw = 20
print("t_pele_froz2tran:".ljust(lw), "%.3E"%t_pele_froz2tran, "[s]")
print("t0:".ljust(lw), "%.3E"%t0, "[s]")
print("t1:".ljust(lw), "%.3E"%t1, "[s]")

print("t_nga_froz:".ljust(lw),  "data_"+"%.3E"%t_nga_froz)
print("1st NGA file:".ljust(lw), re.split("/", fns_sorted[i0])[-1])
print((str(i1-i0+1)+"th NGA file:").ljust(lw), re.split("/", fns_sorted[i1])[-1])

#%%
# Print data
Uin = np.zeros((i1-i0+1,din["ny"],din["nz"]))
Vin = np.zeros((i1-i0+1,din["ny"],din["nz"]))
Win = np.zeros((i1-i0+1,din["ny"],din["nz"]))
if True:
  for ifn, fn in enumerate(fns_sorted[i0:i1+1]):
    print(ifn)
    dataname = fn
    configname = os.path.join(case_path, "config")
    din = NGA_reader.NGA_reader(dataname, configname)
    Uin[ifn,:,:] = din["U"][ix_target,:,:]
    Vin[ifn,:,:] = din["V"][ix_target,:,:]
    Win[ifn,:,:] = din["W"][ix_target,:,:]
np.save("Uin", Uin)
np.save("Vin", Vin)
np.save("Win", Win)

Uin = np.load("Uin.npy")
Vin = np.load("Vin.npy")
Win = np.load("Win.npy")
#%%
nx_out = Uin.shape[0]      # Inlet-outlet (in PeleLMeX)
ny_out = Uin.shape[2] + 1  # Periodic (in PeleLMeX)
nz_out = Uin.shape[1]      # Height (in PeleLMeX)
U = np.zeros((nx_out, ny_out, nz_out))
V = np.zeros((nx_out, ny_out, nz_out))
W = np.zeros((nx_out, ny_out, nz_out))
for j in range(0, ny_out-1):
  for k in range(0, nz_out):
    U[:,j,k] = Uin[:,k,j]
    V[:,j,k] = Vin[:,k,j]
    W[:,j,k] = Win[:,k,j]
# Enforce zero velocity at wall
U[:,:,0] = 0.0; U[:,:,-1] = 0.0
V[:,:,0] = 0.0; V[:,:,-1] = 0.0
W[:,:,0] = 0.0; W[:,:,-1] = 0.0
# Enforce periodic condition along inlet-outlet and trasverse direction
U[:,-1,:] = U[:,0,:]
V[:,-1,:] = V[:,0,:]
W[:,-1,:] = W[:,0,:]
#%%
xf = ts[i0:i1+1] - t_nga_froz
yf = din["z"]
zf = np.zeros((U.shape[2]))
zf[0] = din["y"][1] 
zf[-1] = din["y"][-2]
zf[1:-1] = (din["y"][1:-2]+din["y"][2:-1]) / 2.
zf = zf - zf[0]

folder_out = "%.1E"%t0 + "_" + "%.1E"%t1
if not os.path.exists(folder_out):
  os.makedirs(folder_out)
#%%
Xf, Yf, Zf = np.meshgrid(np.arange(0, xf.shape[0]), 
                         np.arange(0, yf.shape[0]), 
                         np.arange(0, zf.shape[0]), indexing="ij")

Xf1D = Xf.reshape(-1, order="F")
Yf1D = Yf.reshape(-1, order="F")
Zf1D = Zf.reshape(-1, order="F")
U1D = U.reshape(-1, order="F")
V1D = V.reshape(-1, order="F")
W1D = W.reshape(-1, order="F")

data = np.vstack((Xf1D, Yf1D, Zf1D, U1D, V1D, W1D)).T
fmt = "%i, %i, %i, %.18e, %.18e, %.18e"
np.savetxt(os.path.join(folder_out, "channel.uvw"), data, fmt=fmt, delimiter=",", header="i, j, k, u, v, w")
np.savetxt(os.path.join(folder_out, "channel.x"), np.vstack((xf)), header="t")
np.savetxt(os.path.join(folder_out, "channel.y"), np.vstack((yf)), header="y")
np.savetxt(os.path.join(folder_out, "channel.z"), np.vstack((zf)), header="z")

#%%
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
vmin = -10; 
vmax = +10;
cmap = mpl.colormaps["jet"]
im = ax.contourf(xf*10.0, zf, W[:,33,:].transpose(), 
                levels=np.linspace(vmin,vmax,1000), cmap=cmap)
ax.set_xlabel("x", fontsize=20)
ax.set_xlabel("x", fontsize=20)
ax.set_title(r"$W$", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_aspect('equal')
plt.savefig("contourW.png",dpi=300,bbox_inches="tight")
#%%