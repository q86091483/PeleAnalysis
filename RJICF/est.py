#%%
import numpy as np

finest_level = 3
r = np.array([np.mean([1.0]),
              np.mean([0.9]),
              np.mean([0.36]),
              np.mean([0.22])])
tadv = np.mean([170]) / 3600. # [h]
nproc = 2000 # cores
dt_finest = np.mean([9.5E-9])

Lx = 0.0128
Uc = 100.
tf = Lx / Uc

nt_f = tf / dt_finest
nx0 = 256; ny0 = 72; nz0 = 112;
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
#cost_per_cell_step = 3E-7
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