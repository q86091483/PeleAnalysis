#%%
import sys, os
import re
path_PeleAnalysis = os.path.abspath("..")
sys.path.append(path_PeleAnalysis)
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

import cantera as ct
gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/BurkeH2/chem.yaml")

# Input
# Where plt files are stored
plt_folder = "/scratch/b/bsavard/zisen347/scopingRuns/conv3D"
# Case name
case_name = "conv3D"
# Patterns of plotfiles to be processed
plt_pattern = "plt_1*"
# Planes to be extracted
Djet = 5.0E-4
plane_x = np.array([]) * Djet
plane_y = np.array([-2.0]) * Djet
#plane_y = np.array([1.5, 2.0, 3.5]) * Djet
plane_z = np.array([]) * Djet
# Prefix
str_prefix = "HRR_T"
# Max level
max_level = 1
# Fields to be extracted
field_names = ["density", "temp", "mixture_fraction", "mag_vort",
               "x_velocity", "y_velocity", "z_velocity",
               "mixture_fraction_userdef_0", "mixture_fraction_userdef_1",
               "age_0", "age_1"]
field_names = ["density", "temp",
               "x_velocity", "y_velocity", "z_velocity",
               "mixture_fraction_userdef_0", "mixture_fraction_userdef_1",
               "age_0", "age_1", "forcing_age_0", "forcing_age_1",
               "forcing_mixture_fraction_userdef_0", "forcing_mixture_fraction_userdef_1"]
field_names = ["density", "temp","mixture_fraction",
               "x_velocity", "y_velocity", "z_velocity",
               "mixture_fraction_userdef_0", "mixture_fraction_userdef_1",
               "age_0", "age_1", "agepv_0", "agepv_1"]

for isp, spn in enumerate(gas_mix.species_names):
  field_names.append("Y("+spn+")")

# Output data folder
output_dir = "/scratch/b/bsavard/zisen347/PeleAnalysis/Data_age/"
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
output_slice_dir = os.path.join(output_dir, "Slice2D_plt_lev0")
if not os.path.exists(output_slice_dir):
  os.mkdir(output_slice_dir)
output_case_slice_dir = os.path.join(output_slice_dir, case_name)
if not os.path.exists(output_case_slice_dir):
  os.mkdir(output_case_slice_dir)

# Get file names
fnn = os.path.join(plt_folder, plt_pattern)
fns_unsorted = glob.glob(fnn)

def get_key(s):
  ss = re.split("plt_", s)[-1]
  return int(ss)
fns_sorted = sorted(fns_unsorted, key=get_key)
#fns_sorted = fns_unsorted

for fn in fns_sorted:
  print(fn)
#%%

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
