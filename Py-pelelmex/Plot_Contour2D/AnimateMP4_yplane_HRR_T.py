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
# Case name
case_name = "conv3D"
# Patterns of plotfiles to be processed
plt_pattern = "plt_1*"
str_plane = "HRR_T_y=-1.000E-03"
mov_name = "Level1_AGEPV"
zst = 0.0252
Djet = 5.0E-4
# Output data folder
fig_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Figure"
fig_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Figure"
fig_slice_dir = os.path.join(fig_dir, "Slice2D_age_lev0")

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
#%%
if not os.path.exists(fig_dir):
  os.mkdir(fig_dir)
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
#%%
def get_key(s):
  ss = re.split("t=", s)[1]
  sss = ss.split(".png")[0]
  return float(sss)
fns_sorted = sorted(fns_unsorted, key=get_key)
#%%
import cv2
video_name = "./" + mov_name + ".mp4"
frame = cv2.imread(fns_sorted[0])
height, width, layers = frame.shape
#cv2.VideoWriter(output_filename, fourcc, fps, self._window_shape)
video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'XVID'), 10, (width,height))

for image in fns_sorted:
    video.write(cv2.imread(image))


cv2.destroyAllWindows()
video.release()

# %%
