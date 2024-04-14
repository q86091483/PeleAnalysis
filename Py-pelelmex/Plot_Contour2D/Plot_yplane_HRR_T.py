#%%
import sys, os
import re
import numpy as np
import glob
import matplotlib.pyplot as plt

path_PeleAnalysis = os.path.abspath("../..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData

# Input
# Where original plot files are stored
plt_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data/Slice2D/"
              "Micromix/HRR_T_y=9.000E-04") 
# Case name
case_name = "Micromix"
# Patterns of plotfiles to be processed
plt_pattern = "plt_*"
# Field_names
field_names = [[" ", " "], [" ", " "]]
npx = 2; npy = 2
fig_unit_y = 1.5

# Output data folder
fig_dir = "../Figure"
if not os.path.exists(fig_dir):
  os.mkdir(fig_dir)
fig_slice_dir = os.path.join(fig_dir, "Slice2D")
if not os.path.exists(fig_slice_dir):
  os.mkdir(fig_slice_dir)
fig_case_slice_dir = os.path.join(fig_slice_dir, case_name)
if not os.path.exists(fig_case_slice_dir):
  os.mkdir(fig_case_slice_dir)

# Get file names
fns_unsorted = glob.glob(os.path.join(plt_folder, plt_pattern))

def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".npz")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
#%%
for ifn, fn in enumerate(fns_sorted[0:1]):
  output_name = re.split("/", fn)[-1]
  output_name = re.split(".npz", output_name)[0]
  time = float(re.split("t=", output_name)[1])
  output_name = os.path.join(fig_case_slice_dir, output_name)

  f = np.load(fn)
  x1D = f["x"]; Lx = x1D[-1] - x1D[0]
  y1D = f["y"]; Ly = y1D[-1] - y1D[0]
  

  figsize = ((Lx/Ly)*fig_unit_y*npx, fig_unit_y*npx)
  fig, axs = plt.subplots(ncols=npx, nrows=npy, figsize=figsize)
  
  for ipy in range(0, npy):
    for ipx in range(0, npx):
      ax = axs[ipy, ipx]
      if (field_names[ipy][ipx] == "HeatRelease")
        ax
  #print(time)
#%%

