#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
import h5py

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive_debug/Tsave_step593_sdcIter1_deltaTiter9"
ds = yt.load(fn)
dd = ds.all_data()
zst = 0.0252
#%%
xmin = -17.5E-4; xmax = 122.5E-4
ymin = -2E-3; ymax = 2E-3
zmin = -0.0E-3; zmax = 6.0E-3
Lx = xmax - xmin 
Ly = ymax - ymin 
Lz = zmax - zmin

# Normalized parameters
tau_flow = Ly / 20. 

res_x=1000
res_z=res_x * (Lz / Lx)
res_y=res_x * (Ly / Lx)
res_xz = [res_z, res_x]  # create an image with 1000x1000 pixels
res_xy = [res_y, res_x]  # create an image with 1000x1000 pixels
res_yz = [res_y, res_z]  # create an image with 1000x1000 pixels


loc = 2.2188e+00*5E-4
ye = 1.5625e-01 * 5E-4
ze = 1.5938e+00 * 5E-4

fn = "comp0"; dir = "y"; vmin=0; vmax=2000
slc = ds.slice(dir, coord=loc, center=(xmin+0.5*Lx, loc, zmin+0.5*Lz), )
frb = slc.to_frb(width=((Lz, "cm"),(Lx, "cm")), resolution=res_xz)
arr_temp = np.array(frb[fn])
fig, ax = plt.subplots()
im = ax.imshow(arr_temp.transpose(), 
              origin="lower", cmap="viridis", extent=[xmin, xmax, zmin, zmax],
              vmin=vmin, vmax=vmax)
ax.scatter(ye, ze, s=15, color="r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")
#%%