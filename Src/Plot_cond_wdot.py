#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt06524_Y.h5"
f = h5py.File(fn, 'r+')

idx = 3; idy = 0; idz = 0;

nrows = 2; ncols = 2; figunit = 3.5
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(nrows*figunit*1.3,ncols*figunit))
fontsize=22; labelsize=22; ticksize=20; linewidth=2.5

zbin = np.linspace(0.001, 0.13, 101)
spn = "H2O"
for ipx in range(0, nrows):
    for ipy in range(0, ncols):
      ax = axs[ipx, ipy]
      idx = ipx*2 + ipy
      phi = f["DATA"]["uY("+spn+")_mean"][:,idx,idy,idz]
      ax.plot(zbin, phi, linewidth=linewidth, color='b', label="convev")

      #phi_sdr = f["DATA"]["rhoD("+spn+")_mean"][:,idx,idy,idz];
      #ax.plot(zbin, phi_sdr, linewidth=linewidth, color='k', label="SDR") 

      #phi_rr = f["DATA"]["rr("+spn+")_mean"][:,idx,idy,idz];
      #ax.plot(zbin, phi_rr, linewidth=linewidth, color='r', label="SDR") 

      #ax.set_yscale("log")

      ax2 = ax.twinx() 
      phi = f["DATA"]["volume_sum"][:,idx,idy,idz]
      fz = phi / np.sum(phi)
      ax2.bar(zbin, fz, width=zbin[-1]/101, edgecolor=None, color='k', linewidth=linewidth, alpha=0.5)
      ax.set_xlim([zbin[0]-0.004, zbin[-1]+0.004])

      ax.set_title(str(idx), fontsize=fontsize)
      if (ipy == 0):
         #ax.set_yticks([1000, 1500, 2000, 2500])
         s = "s"
         ax.set_ylabel(r"$\langle N_{H2O} | Z \rangle$", fontsize=labelsize)
      #ax.set_yticks([10, 100, 1000])
      
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
