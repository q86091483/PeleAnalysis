#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/Slice2D_derived/Micromix/Derived_y=-1.125E-03/plt_10580_derived_t=6.077E-04.npz"
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3
Djet = 4.5E-04

extent=np.array([xmin, xmax, zmin, zmax]) / Djet
tet = np.load(fn)
hrr = tet["HeatRelease"]
mixfrac = tet["mixture_fraction"]
pv = tet["pv"]
rho = tet["rho"]
FI = tet["FI"]

print("Max/min rho: ",  np.amin(rho), np.amax(rho))
print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min progvar: ",  np.amin(pv), np.amax(pv))
print("Max/min FI: ",  np.amin(FI), np.amax(FI))
#
#%%

if (True):
  # C
  fig, ax = plt.subplots()
  vmin = 0; vmax = 1.2
  im = ax.imshow(pv, origin="lower", 
                 vmin=vmin, vmax=vmax, 
                  cmap="jet", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')

  # hrr
  fig, ax = plt.subplots()
  vmin = -1E11; vmax = 1.0E11
  im = ax.imshow(hrr, origin="lower", 
                  vmin=vmin, vmax=vmax, 
                  cmap="seismic", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')

  # Z 
  fig, ax = plt.subplots()
  #vmin = 0; vmax = 1.0E11
  im = ax.imshow(mixfrac, origin="lower", 
                  #vmin=vmin, vmax=vmax, 
                  cmap="hot", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')

  # FI 
  fig, ax = plt.subplots()
  #vmin = 0; vmax = 1.0E11
  im = ax.imshow(FI, origin="lower", 
                  #vmin=vmin, vmax=vmax, 
                  cmap="seismic", extent=extent)
  divider = make_axes_locatable(ax)
  #ax.contour(pv, levels=[0.0252],
  #          origin='lower', 
  #          colors=['white'], extent=extent)
  cax = divider.append_axes('right', size='5%', pad=0.05)
  fig.colorbar(im, cax=cax, orientation='vertical')
  #%%




#%%