#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)
import matplotlib.pyplot as plt
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/Slice2D/Micromix/Derived_y=9.000E-04/plt_07000_derived_t=0.000E+00.npz"
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3

extent=np.array([xmin, xmax, zmin, zmax]) / Djet
tet = np.load(fn)
hrr = tet["HeatRelease"]
mixfrac = tet["mixture_fraction"]
pv = tet["pv"]


print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min prog ar: ",  np.amin(pv), np.amax(pv))



if (True):
  fig, ax = plt.subplots()
  vmin = 0; vmax = 1.0
  ax.imshow(pv, vmin=vmin, vmax=vmax, origin="lower", 
            cmap="jet", extent=extent)
  ax.contour(mixfrac, levels=[0.0252],
            origin='lower', 
            colors=['white'], extent=extent)


#%%