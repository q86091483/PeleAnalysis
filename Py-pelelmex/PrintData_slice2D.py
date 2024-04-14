#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob

# Input
# Where plt files are stored
plt_folder = "/scratch/b/bsavard/zisen347/scopingRuns/MicroMix"
# Case name
case_name = "Micromix"
# Patterns of plotfiles to be processed
plt_pattern = "plt_0[4,5,6]*"
# Planes to be extracted
plane_x = []
plane_y = [-9E-4, 9E-4]
plane_z = [2E-3, 3E-3]
# Prefix
str_prefix = "HRR_T"
# Fields to be extracted
field_names = ["HeatRelease", "temp", "mixture_fraction", "Y(NO)", "Y(N2O)", "Y(NNH)"]
# Max level
max_level = 1
# Output data folder
output_dir = "./Data"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
output_slice_dir = os.path.join(output_dir, "Slice2D")
if not os.path.exists(output_slice_dir):
  os.mkdir(output_slice_dir)
output_case_slice_dir = os.path.join(output_slice_dir, case_name)
if not os.path.exists(output_case_slice_dir):
  os.mkdir(output_case_slice_dir)

# Get file names
fns_unsorted = glob.glob(os.path.join(plt_folder, plt_pattern))

def get_key(s):
  ss = re.split("plt_", s)[-1]
  return int(ss)
fns_sorted = sorted(fns_unsorted, key=get_key)

for iy, pos in enumerate(plane_y):
  normal = 1

  folder_name = str_prefix + "_y=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):   
    os.mkdir(folder_name)
  
  for ifn, fn in enumerate(fns_sorted):
    out_name = re.split("/", fn)[-1]
    time = HeaderData(fn).time
    out_name = out_name + "_t=" + "%.3E"%time
    out_name = os.path.join(folder_name, out_name)

    mand = Mandoline(fn, 
                     fields=field_names, 
                     limit_level=max_level,
                     serial=True,
                     verbose=1)
    mand.slice(normal=normal, 
               pos=pos,
               outfile=out_name,
               fformat="array",
               uselog=1,
               )

for iz, pos in enumerate(plane_z):
  normal = 1

  folder_name = str_prefix + "_z=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):   
    os.mkdir(folder_name)
  
  for ifn, fn in enumerate(fns_sorted):
    out_name = re.split("/", fn)[-1]
    time = HeaderData(fn).time
    out_name = out_name + "_t=" + "%.3E"%time
    out_name = os.path.join(folder_name, out_name)

    mand = Mandoline(fn, 
                     fields=field_names, 
                     limit_level=max_level,
                     serial=True,
                     verbose=1)
    mand.slice(normal=normal, 
               pos=pos,
               outfile=out_name,
               fformat="array",
               uselog=1,
               )
#%%
