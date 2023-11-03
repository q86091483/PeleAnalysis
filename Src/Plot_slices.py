#%%
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import yt

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Load data
plt_name = "4840"; time = 000
time_name = "$t="+str(time)+"\;[\mathrm{\mu s}]$"
fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt0"+plt_name+"_T_Z"
ds = yt.load(fn)
dd = ds.all_data()

fn = "/scratch/b/bsavard/zisen347/PeleAnalysis/Src/cond_ISRN/plt0"+plt_name+"_Y.h5"
f = h5py.File(fn, 'r+')
nzbin = 101
zbin = np.linspace(0.001, 0.13, nzbin)

Lx = 0.0033; Ly = 0.0009; Lz = 0.0012
res_x=500
res_y=res_x * (Ly / Lx)
res = [res_x, res_y]  # create an image with 1000x1000 pixels

slc = ds.slice("z", coord=Lz/4, center=(Lx/2, Ly/2, Lz/4), )
frb = slc.to_frb(width=((Lx, "cm"),(Ly, "cm")), resolution=res)
#frb = proj.to_frb(width=(0.00165, "cm"), resolution=res, center=(0.00165, 0.00045, 0.0003), )
arr_temp = np.array(frb["temp"])

labelsize = 20; ticksize = 18; titlesize = 20; linewidth = 2.5;
figg = plt.figure(constrained_layout=True)
figunit = 1; ncols = 3; nrows = 2;
widths = [figunit * 1.2] * ncols
heights = [figunit] * nrows
spec = figg.add_gridspec(ncols=3, nrows=2, width_ratios=widths,
                         height_ratios=heights)

irow = 0
ax = figg.add_subplot(spec[0, :])
ax.imshow(arr_temp, origin="lower", cmap="turbo", vmax=2700, vmin=850, extent=[0,Lx,0,Ly])
ax.set_xticks([0, 1E-3, 2E-3, 3E-3])
ax.set_xticklabels(["0.0", "1.0", "2.0", "3.0"])
ax.set_xlabel("$x$ [cm]", fontsize=labelsize)
ax.xaxis.set_label_coords(.5, -.1)
ax.set_yticks([0, 0.45E-3, 0.9E-3])
ax.set_yticklabels(["0.0", "0.45", "0.9"])
ax.set_ylabel("$y$ [mm]", fontsize=labelsize)
ax.tick_params(axis='both', which='major', labelsize=ticksize)
ax.tick_params(axis='both', which='minor', labelsize=ticksize)
ax.set_title(time_name, fontsize=titlesize)
ax.plot([1.1E-3, 1.1E-3], [0.0, 0.9E-3], linestyle="--", linewidth=3.5, color="w")
ax.plot([2.2E-3, 2.2E-3], [0.0, 0.9E-3], linestyle="--", linewidth=3.5, color="w")

spn = "H2O"
irow = 1; idy = 0; idz = 0
for icol in range(len(widths)):
   ax1 = figg.add_subplot(spec[irow, icol])
   idx = icol
   phi = f["DATA"]["temp_mean"][:,idx,idy,idz]
   ax1.plot(zbin, phi/2700, linewidth=linewidth, color='r')
   ax1.set_ylim([-1.0, 1.0])

   ax2 = ax1.twinx() 
   phi = f["DATA"]["volume_sum"][:,idx,idy,idz]
   fz = phi / np.sum(phi)
   ax2.bar(zbin, fz, width=zbin[-1]/nzbin, edgecolor=None, color='k', linewidth=linewidth, alpha=0.5)
   ax2.set_ylim([0, 0.06])
   if (idx == 2):
      ax2.set_ylabel(r"$f_{Z,\Delta}$", fontsize=labelsize)
      ax2.set_yticks([0, 0.06])
      ax2.yaxis.set_label_coords(1.08, 0.5)
   else:
      ax2.set_yticklabels([])

   phi_sdr = f["DATA"]["rhoD("+spn+")_mean"][:,idx,idy,idz]
   phi_sdr = phi_sdr * f["DATA"]["Y("+spn+")_mean"][:,idx,idy,idz]
   ax1.plot(zbin, phi_sdr / np.amax(phi_sdr), linewidth=linewidth, color='k', label="SDR")

   phi_conv = f["DATA"]["uY("+spn+")_mean"][:,idx,idy,idz]
   ax1.plot(zbin, phi_conv / np.amax(phi_conv), linewidth=linewidth, color='b', label="convev")

   if (idx == 0):
      ax1.set_ylabel(r"$\langle \phi | Z \rangle_\Delta/\langle \phi | Z \rangle_\Delta,max$", fontsize=labelsize)
   ax1.tick_params(axis='both', which='major', labelsize=ticksize)
   ax1.tick_params(axis='both', which='minor', labelsize=ticksize)
   ax2.tick_params(axis='both', which='major', labelsize=ticksize-2)
   ax2.tick_params(axis='both', which='minor', labelsize=ticksize-2)
   ax1.set_xlabel("$Z$", fontsize=labelsize)

plt.savefig(str(time)+".png", dpi=300, bbox_inches="tight")



#%%