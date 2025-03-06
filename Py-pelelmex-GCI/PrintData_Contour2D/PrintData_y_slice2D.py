#%%
import sys, os
import re
path_PeleAnalysis = "/scratch/w47/zl5403/PeleAnalysis"
sys.path.append(path_PeleAnalysis)

from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

import cantera as ct

gas_mix = ct.Solution(
            os.path.join(path_PeleAnalysis,
            "Py-pelelmex-GCI/Input/prf_rd53/qssa.yaml"))

# Input
# Where plt files are stored
case_name = "Ujet30_Re7200_rd53"
#case_name = "GCI_jet3D"
#plt_folder = "/scratch/b/bsavard/zisen347/scopingRuns/NUIG_Re4000_2J6_4atm/Level_3"
plt_folder = "/scratch/w47/zl5403/scopingRuns/" + case_name

# Case name
# Patterns of plotfiles to be processed
plt_pattern = "plt_07125*"
# Planes to be extracted
Djet = 1.7E-4
plane_x = np.array([]) * Djet
plane_y = np.array([0.0]) * Djet
plane_z = np.array([]) * Djet
# Prefix
str_prefix = "HRR_T"
# Max level
max_level = 1
# Fields to be extracted
field_names = ["density", "temp", "mag_vort",
               "x_velocity", "y_velocity", "z_velocity", "mixture_fraction",
               "mixture_fraction_userdef_0", "mixture_fraction_userdef_1",
               "age_0", "age_1"]
#field_names = ["density", "temp", "mixture_fraction", "mag_vort",
#               "x_velocity", "y_velocity", "z_velocity",]
#               ["age_0", "age_1", "mixture_fraction_userdef_0", "mixture_fraction_userdef_1"]]
if True:
  species_names = ['N2','H2','H','O2','O','H2O','OH','H2O2',
                   'HO2','CO','CO2','CH4','CH3','CH2','CH2O','HCO',
                   'HCOH','C2H6', 'C2H5','C2H5O2H', 'C2H5O2','C2H4','C2H3','CHOCHO',
                   'C2H3OO','CHCHO','C2H2','PC2H4OH', 'O2C2H4OH','CH2CHO','CH2CO','HCCO',
                   'C3H8','IC3H7','C3H6','C3H5-A','CH3CHCHO','AC3H5OOH','C3H4-P','C3H4-A',
                   'C3H3','C3H2','C2H3CHO','NC7H16', 'C7H15O2','C7KET','C5H11CO','IC8H18', 'C8H17O2','C8KET','C6H13CO','C8H16',
                   'C6H6']
  for isp, spn in enumerate(species_names):
    field_names.append("Y("+spn+")")
    print(isp, spn)
else:
  species_names = gas_mix.species_names
  for isp, spn in enumerate(gas_mix.species_names):
    field_names.append("Y("+spn+")")
    print(isp, spn)
#%%
# Output data folder
output_dir = os.path.join(path_PeleAnalysis, "Data")
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
output_slice_dir = os.path.join(output_dir, "Slice2D_plt_lev0")
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
