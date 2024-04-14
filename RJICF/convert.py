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

# Inputs
case_folder   = "/scratch/b/bsavard/zisen347/cases/"
case_name     = "NUIG_Re4000_2J6_4atm"
fields = ["P"] * 6 + ["V"] * 6 + ["RHO"] * 6
idirs = [1, 1, 2, 2, 3, 3] * 3
isls = [1, 96, 1, 96, 1, 96] * 3
case_path     = os.path.join(case_folder, case_name)

dataname = os.path.join(case_path, "data_1.000E-02")
#dataname = os.path.join(case_path, "data.init")
configname = os.path.join(case_path, "config")
din = NGA_reader.NGA_reader(dataname, configname)

nx_out = din["U"].shape[0] + 1  # Inlet-outlet (in PeleLMeX)
ny_out = din["U"].shape[2] + 1  # Periodic (in PeleLMeX)
nz_out = din["U"].shape[1]      # Height (in PeleLMeX)
U = np.zeros((nx_out, ny_out, nz_out))
V = np.zeros((nx_out, ny_out, nz_out))
W = np.zeros((nx_out, ny_out, nz_out))
for j in range(0, ny_out-1):
  for k in range(0, nz_out):
    U[0:-1,j,k] = din["U"][:,k,j]
    V[0:-1,j,k] = din["V"][:,k,j]
    W[0:-1,j,k] = din["W"][:,k,j]
# Enforce zero velocity at wall
U[:,:,0] = 0.0; U[:,:,-1] = 0.0
V[:,:,0] = 0.0; V[:,:,-1] = 0.0
W[:,:,0] = 0.0; W[:,:,-1] = 0.0
# Enforce periodic condition along inlet-outlet and trasverse direction
U[:,-1,:] = U[:,0,:]
V[:,-1,:] = V[:,0,:]
W[:,-1,:] = W[:,0,:]
U[-1,:,:] = U[0,:,:]
V[-1,:,:] = V[0,:,:]
W[-1,:,:] = W[0,:,:]

xf = din["x"]
yf = din["z"]
zf = np.zeros((U.shape[2]))
zf[0] = din["y"][1] 
zf[-1] = din["y"][-2]
zf[1:-1] = (din["y"][1:-2]+din["y"][2:-1]) / 2.
zf = zf - zf[0]

if (False):
  for j in range(0, yf.shape[0]):
    for k in range(0, zf.shape[0]):
      #U[:,j,k] = k
      V[:,j,k] = k/100
      #W[:,j,k] = 0.0

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
np.savetxt("channel.uvw", data, fmt=fmt, delimiter=",", header="i, j, k, u, v, w")
np.savetxt("channel.x", np.vstack((xf)), header="x")
np.savetxt("channel.y", np.vstack((yf)), header="y")
np.savetxt("channel.z", np.vstack((zf)), header="z")

#%%
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
vmin = -7; 
vmax = 7;
cmap = mpl.colormaps["jet"]
im = ax.contourf(xf, zf, W[:,33,:].transpose(), 
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
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
vmin = -7; 
vmax = 7;
cmap = mpl.colormaps["jet"]
im = ax.contourf(xf, zf, V[:,33,:].transpose(), 
                levels=np.linspace(vmin,vmax,1000), cmap=cmap)
ax.set_xlabel("x", fontsize=20)
ax.set_xlabel("x", fontsize=20)
ax.set_title(r"$V$", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_aspect('equal')
plt.savefig("contourV.png",dpi=300,bbox_inches="tight")

# %%
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
vmin = 0.; 
vmax = 45;
cmap = mpl.colormaps["jet"]
im = ax.contourf(xf, zf, U[:,33,:].transpose(), 
                levels=np.linspace(vmin,vmax,1000), cmap=cmap)
ax.set_xlabel("x", fontsize=20)
ax.set_xlabel("x", fontsize=20)
ax.set_title(r"$U $", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_aspect('equal')
plt.savefig("contourU.png",dpi=300,bbox_inches="tight")

# %%
