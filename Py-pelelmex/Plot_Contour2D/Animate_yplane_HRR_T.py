#%%
import sys, os
import re
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt

path_PeleAnalysis = os.path.abspath("../..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData

# Input
# Where original plot files are stored
fig_folder = ("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/"
              "Figure/Slice2D/Micromix")
# Case name
case_name = "Micromix"
# Patterns of plotfiles to be processed
plt_pattern = "plt_*"
str_plane = "HRR_T_y=9.000E-04"
# Field_names
field_names = [["HeatRelease", "temp"], ["Y(NO)", "Y(N2O)"]]
# 
zst = 0.0252
Djet = 4.5E-4
# Plot parameter
npx = 2; npy = 2
fig_unit_y = 1.5
labelsize = 16
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
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
fig_plane_case_slice_dir = os.path.join(fig_case_slice_dir, str_plane)
if not os.path.exists(fig_plane_case_slice_dir):
  os.mkdir(fig_plane_case_slice_dir)

# Get file names
fns_unsorted = glob.glob(os.path.join(fig_plane_case_slice_dir, plt_pattern))

def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".png")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
#%%

import imageio
images = []
for filename in fns_sorted:
    print(filename)
    images.append(imageio.imread(filename))
imageio.mimsave("./" + str_plane+ "_movie.gif", images)
#%%