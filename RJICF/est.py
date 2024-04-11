#%%
import numpy as np

finest_level = 1
r = np.array([np.mean([1.0]), 
              np.mean([1.0]), 
              np.mean([0.7916]), 
              np.mean([0.4011, 0.4022, 0.4059, 0.4055, 0.4153 ]), 
              np.mean([0.05])  ])
r = np.array([np.mean([1.0]),
              np.mean([0.833]),
              np.mean([0.4]),
              np.mean([0.05])])
tadv = np.mean([82.83]) / 3600. # [h]
nproc = 160 # cores
dt_finest = np.mean([5.688723736e-08*2])

Lx = 0.014
Uc = 35.583
tf = Lx / Uc

nt_f = tf / dt_finest
nx0 = 336; ny0 = 96; nz0 = 144;
n0 = nx0 * ny0 * nz0
max_level = r.shape[0]

def get_n(lev, r):
  n = nx0 * ny0 * nz0 * r * np.power(8, lev)
  return n

nlev = []
print("Grid # at level:")
for i in range(0, max_level):
  nlev.append(get_n(i, r[i]))
  print(i, ": ", nlev[i])

total_cell_finest = 0; 
for lev in range(0, finest_level+1):
  total_cell_finest = total_cell_finest + nlev[lev]

cost_per_cell_step = nproc * tadv / total_cell_finest
cost_per_cell_step = 3E-7
print("cost_per_cell_step", cost_per_cell_step, "\n")

nlast = 0
for i, lev in enumerate(range(0, finest_level)):
  nlast = nlast + get_n(lev, r[lev])

for i, lev in enumerate(range(finest_level, max_level)):
  Ni = (tf / dt_finest) * np.power(2, lev - finest_level)
  ni = get_n(lev, r[lev])
  print("Level ", lev)
  print("Ni", Ni) 
  print("ni: ", ni)

  nlast = ni + nlast
  cost = Ni * nlast * cost_per_cell_step
  print("nlast", nlast)
  print("Cost lev", cost, "\n")
#%%