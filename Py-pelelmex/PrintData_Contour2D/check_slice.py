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

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/Slice2D/Micromix/HRR_T_y=1.125E-03/plt_10080_t=5.998E-04.npz"
tet = np.load(fn)
x_velocity = tet["temp"]
y_velocity = tet["y_velocity"]
z_velocity = tet["z_velocity"]

print("Max/min x_velocity: ",  np.amin(x_velocity), np.amax(x_velocity))
print("Max/min y_velocity: ",  np.amin(y_velocity), np.amax(y_velocity))
print("Max/min z_velocity: ",  np.amin(z_velocity), np.amax(z_velocity))

if (True):
  fig, ax = plt.subplots()
  ax.imshow(x_velocity, vmin=-100, vmax=2800, origin="lower")

  fig, ax1 = plt.subplots()
  ax1.imshow(y_velocity, vmin=-100, vmax=100, origin="lower")

  fig, ax2 = plt.subplots()
  ax2.imshow(z_velocity, vmin=-100, vmax=100, origin="lower")

#%%