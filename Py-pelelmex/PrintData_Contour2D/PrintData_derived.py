#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

# Input
# Where plt files are stored
plt_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/res_MicroMix_RPA/"
# Case name
case_name = "Micromix"
# Patterns of plotfiles to be processed
plt_pattern = "plt_10580_derived"
# Planes to be extracted
Djet = 4.5E-4
#plane_x = np.array([1.0, 2.0, 5.0, 8.0, 15.0, 24.0]) * Djet
plane_x = np.array([]) * Djet
#plane_y = np.array([-3.5, -2.5, -1.5, 1.5, 2.0, 2.5, 3.5]) * Djet
plane_y = np.array([-2.0, -2.5]) * Djet
plane_z = np.array([]) * Djet
# Prefix
str_prefix = "Derived"
# Fields to be extracted
# Max level
max_level = 1
field_names = ["rho", "HeatRelease", "mixture_fraction", "temp", "Y(H2)", "pv", "FI", 
               "HeatReleaseFI", "rhorr(NO)", "rhorr(NNH)", "rhorr(N2O)",
               "R10", "zone"]
#field_names = ["x_velocity", "y_velocity", "z_velocity"]
# Output data folder
output_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Data"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
output_slice_dir = os.path.join(output_dir, "Slice2D_derived")
if not os.path.exists(output_slice_dir):
  os.mkdir(output_slice_dir)
output_case_slice_dir = os.path.join(output_slice_dir, case_name)
if not os.path.exists(output_case_slice_dir):
  os.mkdir(output_case_slice_dir)

# Get file names
fns_unsorted = glob.glob(os.path.join(plt_folder, plt_pattern))

def get_key(s):
  ss = re.split("_derived", s)[0]
  ss = re.split("plt_", ss)[-1]
  return int(ss)
fns_sorted = sorted(fns_unsorted, key=get_key)
#%%
for fn in fns_sorted:
  print(fn)
for ix, pos in enumerate(plane_x):
  normal = 0

  folder_name = str_prefix + "_x=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):   
    os.mkdir(folder_name)
  
  for ifn, fn in enumerate(fns_sorted):
    print("Processing ix =", ix, ", pos=", pos)

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
for iy, pos in enumerate(plane_y):
  normal = 1

  folder_name = str_prefix + "_y=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):   
    os.mkdir(folder_name)
  
  for ifn, fn in enumerate(fns_sorted):
    print("Processing iy =", iy, ", pos=", pos)

    normal = 1
    out_name = re.split("/", fn)[-1]
    time = HeaderData(fn).time
    out_name = out_name + "_t=" + "%.3E"%time
    out_name = os.path.join(folder_name, out_name)

    print(normal, fn, out_name)
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
for iz, pos in enumerate(plane_z):
  normal = 2

  folder_name = str_prefix + "_z=" + "%.3E" % pos
  folder_name = os.path.join(output_case_slice_dir, folder_name)
  if not os.path.exists(folder_name):   
    os.mkdir(folder_name)
  
  for ifn, fn in enumerate(fns_sorted):
    print("Processing iz =", iz, ", pos=", pos)

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
