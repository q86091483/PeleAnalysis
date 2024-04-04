#%%
import numpy as np


vol_ratios = np.array([[1.0, 0.773, 0.35],
                ])
finest_dts = np.array([5.787635397e-08])

finest_level = 0;
max_level = 3;
finest_dt = np.mean(finest_dts);

num_cells = np.zeros((max_level)) - 1E10; 

nx0 = 224;
ny0 = 64;
nz0 = 96;

for lev in range(0, finest_level+1):
  nx = np.power(nx0, lev)
  ny = np.power(ny0, lev)
  nz = np.power(nz0, lev)
  r = np.mean(vol_ratios, axis=0)[lev]
  num_cells[lev] = nx0 * ny0 * nz0 * np.power(np.power(2, lev),3) * r

nproc = 160
tf = 0.014/17.5
cost_dt = 31 / 3600.
cost_per_second_cell = nproc * cost_dt / np.sum(num_cells) /finest_dt
print("CPU hours per second per cell: ", cost_per_second_cell)
print("CPU hours total", cost_per_second_cell*tf*np.sum(num_cells) )

#with 8519680 cells,2080 boxes, over 77.38095238% of the domain 
#Making new level 2 from scratch
#with 18874880 cells,6702 boxes, over 21.42915272% of the domain 
#%%