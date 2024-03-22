#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt_11000_Y.h5"
#fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt_11000_Y.h5"
f = h5py.File(fn, 'r+')

outNames = ["rho", "covered", "ts11", "ts22", "ts33", "ts12", "ts13", "ts23", "ts"]


wt = np.array(f["DATA"]["volume_mean"])

rho = np.array(f["DATA"]["rho_mean"]) / wt

ts11 = np.array(f["DATA"]["ts11_mean"]) / wt
ts22 = np.array(f["DATA"]["ts22_mean"]) / wt
ts33 = np.array(f["DATA"]["ts33_mean"]) / wt
ts12 = np.array(f["DATA"]["ts12_mean"]) / wt
ts13 = np.array(f["DATA"]["ts13_mean"]) / wt
ts23 = np.array(f["DATA"]["ts23_mean"]) / wt
mu = np.array(f["DATA"]["mu_mean"]) / wt

ts_sum = np.array(f["DATA"]["ts_sum_mean"]) / wt
sum_of_meanprod= ts11*ts11 + ts22*ts22 + ts33*ts33+ \
                 2*ts12*ts12 + 2*ts13*ts13 + 2*ts23*ts23
eps = ts_sum - sum_of_meanprod
eps = 2 * eps * mu / rho

eta = np.power(mu/rho, 3) / eps
eta = np.power(eta, 0.25)

fig, ax = plt.subplots()
vmin = 1E-7; vmax = 1E2
im = ax.imshow(eta[:,5,:].transpose(), cmap="tab20",
               origin="lower",
               norm=LogNorm(vmin, vmax))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title("${\eta_k} \; [\mathrm{m}]$", fontsize=20)
#plt.savefig(fn+"_"+dir+".png", dpi=300, bbox_inches="tight")
#phi = sum_varwt[:,:,:u / sum_wt[:,:,:]
#%%
idx = 3; idy = 0; idz = 0;

nrows = 2; ncols = 2; figunit = 3.5
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(nrows*figunit*1.3,ncols*figunit))
fontsize=22; labelsize=22; ticksize=20; linewidth=2.5

zbin = np.linspace(0.001, 0.13, 101)
for ipx in range(0, nrows):
    for ipy in range(0, ncols):
      ax = axs[ipx, ipy]
      idx = ipx*2 + ipy
      phi = f["DATA"]["uY(H2O)_mean"][:,idx,idy,idz]
      ax.plot(zbin, phi, linewidth=linewidth, color='r')

      ax2 = ax.twinx() 
      phi = f["DATA"]["volume_sum"][:,idx,idy,idz]
      fz = phi / np.sum(phi)
      ax2.bar(zbin, fz, width=zbin[-1]/101, edgecolor=None, color='k', linewidth=linewidth, alpha=0.5)
      ax.set_xlim([zbin[0]-0.004, zbin[-1]+0.004])

      ax.set_title(str(idx), fontsize=fontsize)
      if (ipy == 0):
         #ax.set_yticks([1000, 1500, 2000, 2500])
         ax.set_ylabel(r"$\langle T| Z \rangle ~[K]$", fontsize=labelsize)
      else:
         ax.set_yticks([])
      
      ax2.set_ylim([0., np.amax(fz)])
      if (ipy == ncols-1):
         ax2.set_ylabel(r"$P(Z)$", fontsize=labelsize)  
      if (ipx == 1):
        ax.set_xlabel(r"$Z$", fontsize=labelsize)
      ax.set_xticks([0.0, 0.05, 0.10, 0.13])
      #ax.set_ylim([900, 2700])

      ax.tick_params(axis='both', which='major', labelsize=ticksize)
      ax.tick_params(axis='both', which='minor', labelsize=ticksize)
      ax2.tick_params(axis='both', which='major', labelsize=ticksize)
      ax2.tick_params(axis='both', which='minor', labelsize=ticksize)
    
fig.tight_layout(pad=1.5)
plt.savefig("T_Z.png", dpi=300, bbox_inches="tight")
#ax2 = ax.twinx()
#volume_sum = f["DATA"]["volume_sum"][:,idx,idy,idz]
#ax2.plot(volume_sum / np.sum(volume_sum))

# %%
