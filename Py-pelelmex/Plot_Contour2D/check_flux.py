#%%
import sys, os
import re
path_PeleAnalysis = "/scratch/b/bsavard/zisen347/PeleAnalysis"
sys.path.append(path_PeleAnalysis)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from amr_kitchen.mandoline import Mandoline
from amr_kitchen import HeaderData
import glob
import numpy as np

import cantera as ct
#gas_mix = ct.Solution("/scratch/b/bsavard/zisen347/PeleAnalysis/Py-pelelmex/Input/nuig_H2_4atm/chem.yaml")
gas_mix = ct.Solution("/home/b/bsavard/zisen347/codes/exascale/PelePhysics/Support/Mechanism/Models/nuig_H2_4atm_CHNO/nuig_H2_32sp_4atm.yaml", transport_model="Multi")
Nsp = len(gas_mix.species_names)
H2_ID = gas_mix.species_index("H2")
O2_ID = gas_mix.species_index("O2")
N2_ID = gas_mix.species_index("N2")
NO_ID = gas_mix.species_index("NO")
N2O_ID = gas_mix.species_index("N2O")
NNH_ID = gas_mix.species_index("NNH")

data_folder = "/scratch/b/bsavard/zisen347/PeleAnalysis/Data/"
fn_derived = os.path.join(data_folder, "Slice2D_derived_lev1/MicroMix_age/Derived_y=9.000E-04/plt_17300_derived_t=9.998E-04.npz")
fn_plt     = os.path.join(data_folder, "Slice2D_plt_lev1/MicroMix_age/HRR_T_y=9.000E-04/plt_17300_t=9.998E-04.npz")
xmin = -15.75E-4; xmax = 112.25E-4
ymin = -1.8E-3; ymax = +1.8E-3
zmin = 0.0; zmax = 5.6E-3
Djet = 4.5E-04
extent=np.array([xmin, xmax, zmin, zmax]) / Djet
#%%
# Read data
P0      = 101325*4.0
tet_plt = np.load(fn_plt)
tet_der = np.load(fn_derived)
temp    = tet_der["temp"]
density = tet_plt["density"]
velx    = tet_plt["x_velocity"]
vely    = tet_plt["y_velocity"]
velz    = tet_plt["z_velocity"]
pres    = np.zeros_like(density) + P0
hrr     = tet_plt["HeatRelease"]
hrrFI   = tet_der["HeatReleaseFI"]
mixfrac = tet_der["mixture_fraction"]
pv      = tet_der["pv"]
rho     = tet_der["rho"]
rhorr_NO = tet_der["rhorr(NO)"]
rhorr_N2O = tet_der["rhorr(N2O)"]
rhorr_NNH = tet_der["rhorr(NNH)"]


FI      = tet_der["FI"]
R10     = tet_der["R10"]
T10     = tet_der["T10"]
D10     = tet_der["D10"]
W10     = tet_der["W10"]
C10     = tet_der["C10"]

print("Max/min rho: ",  np.amin(rho), np.amax(rho))
print("Max/min HeatRelease: ",  np.amin(hrr), np.amax(hrr))
print("Max/min Mixfrac: ",  np.amin(mixfrac), np.amax(mixfrac))
print("Max/min progvar: ",  np.amin(pv), np.amax(pv))
print("Max/min FI: ",  np.amin(FI), np.amax(FI))
nz = temp.shape[0]
nx = temp.shape[1]

# Array of species
tet_Y = np.zeros((nz, nx, Nsp))
for isp in range(0, Nsp):
  spn = gas_mix.species_names[isp]
  tet_Y[:,:,isp] = tet_plt["Y("+spn+")"][:,:]

# Load data into ct.SolutionArray
fstate = ct.SolutionArray(gas_mix, (nz, nx))
for i in range(0, nz):
  fstate[i].TPY = temp[i,:], pres[i,:], tet_Y[i,:,:]

# Way 1
net_rates_of_progress = fstate.net_rates_of_progress
forward_rates_of_progress = fstate.forward_rates_of_progress
reverse_rates_of_progress = fstate.reverse_rates_of_progress

net_rates_of_progress_ddC = fstate.net_rates_of_progress_ddC
delta_enthalpy = fstate.delta_enthalpy
hrr_state = np.sum(net_rates_of_progress * delta_enthalpy, axis=(2))
# Way 2
net_production_rates = fstate.net_production_rates
partial_molar_enthalpies = fstate.partial_molar_enthalpies
hrr_state = np.sum(net_production_rates * partial_molar_enthalpies, axis=(2))

#%%
# Reaction fstate
ir = 71
vmin = -1E-2
vmax = 1E-2
fig, ax = plt.subplots()
im = ax.imshow(R10[:,:]/1000., origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$Pele: " + str(gas_mix.reactions()[ir]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower',
#          colors=['white'], ext0ent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
im = ax.imshow(fstate.net_rates_of_progress[:,:,ir], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$Cantera: " + str(gas_mix.reactions()[ir]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower',
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
im = ax.imshow(R10[:,:]/1000-fstate.net_rates_of_progress[:,:,ir], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$Pele-Cantera: " + str(gas_mix.reactions()[ir]) + " [kmol/m^3]$")
divider = make_axes_locatable(ax)
#ax.contour(pv, levels=[0.0252],
#          origin='lower',
#          colors=['white'], extent=extent)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')


#%% Diffusivity coefficient
fig, ax = plt.subplots()
vmin = 0; vmax = 2E-2
im = ax.imshow(T10[:,:], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="jet", extent=extent)
ax.set_title(r"$\rho D_\mathrm{NNH}~\mathrm{[kg/m \cdot s]}" + "-Pele $", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
isp = gas_mix.species_index("NNH")
im = ax.imshow(fstate.mix_diff_coeffs[:,:,isp]*fstate.density, origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="jet", extent=extent)
ax.set_title(r"$\rho D_\mathrm{NNH}~\mathrm{[kg/m\cdot s]}" + "-Cantera$", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

def grad_x(phi):
  dx = (xmax - xmin) / (nx - 1)
  res = np.gradient(phi, axis = (1)) / dx
  return res
def grad_z(phi):
  dz = (zmax - zmin) / (nz - 1)
  res = np.gradient(phi, axis = (0)) / dz
  return res

#%% Diffusion component
fig, ax = plt.subplots()
vmin = -5E1; vmax = -vmin
im = ax.imshow(D10[:,:], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$Y_\mathrm{H:}~\nabla_x \cdot (\rho D \nabla_x ) + \nabla_z \cdot (\rho D \nabla_z)- Pele$", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

isp = gas_mix.species_index("H")
phi = tet_plt["Y(H)"]
gradx_Y = grad_x(phi)
gradz_Y = grad_z(phi)
rhoDgradx = fstate.mix_diff_coeffs[:,:,isp]*fstate.density*gradx_Y
rhoDgradz = fstate.mix_diff_coeffs[:,:,isp]*fstate.density*gradz_Y
lapx = grad_x(rhoDgradx)
lapz = grad_z(rhoDgradz)
fig, ax = plt.subplots()
im = ax.imshow(lapx + lapz, origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$Y_\mathrm{H:}~\nabla_x \cdot (\rho D \nabla_x ) + \nabla_z \cdot (\rho D \nabla_z)- Cantera$", fontsize = 18)

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

#%% convection
fig, ax = plt.subplots()
vmin = -5E1; vmax = -vmin
im = ax.imshow(C10[:,:], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$ Y_\mathrm{H}:~\nabla_x \cdot ( \rho u_x Y ) + \nabla_z \cdot ( \rho u_z Y ) - Pele$",
             fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

phi = tet_plt["Y(H)"]
convY_x = grad_x(density * velx * phi)
convY_z = grad_z(density * velz * phi)
fig, ax = plt.subplots()
isp = gas_mix.species_index("H")
im = ax.imshow(convY_x + convY_z, origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$ Y_\mathrm{H}:~\nabla_x \cdot ( \rho u_x Y ) + \nabla_z \cdot ( \rho u_z Y ) - Cantera$",
             fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

#%% DW
fig, ax = plt.subplots()
vmin = -5E1; vmax = -vmin
im = ax.imshow(W10[:,:], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$ Y_\mathrm{H}:~\nabla_x \cdot ( \rho D \frac{\nabla_z \overline{W}}{\overline{W}}) + \nabla_z \cdot ( \rho D \frac{\nabla_z \overline{W}}{\overline{W}}) "
              "- Pele $", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

fig, ax = plt.subplots()
isp = gas_mix.species_index("H")
Wbar = fstate.mean_molecular_weight / 1000
gradW_x = grad_x(Wbar)
gradW_z = grad_z(Wbar)
Yt = tet_plt["Y(H)"]
phi = Yt * fstate.mix_diff_coeffs[:,:,isp]*fstate.density * gradW_x / Wbar
lapW_x = grad_x(phi)
phi = Yt * fstate.mix_diff_coeffs[:,:,isp]*fstate.density * gradW_z / Wbar
lapW_z = grad_z(phi)
im = ax.imshow(lapW_x + lapW_z, origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$ Y_\mathrm{H}:~\nabla_x \cdot ( \rho D \frac{\nabla_z \overline{W}}{\overline{W}}) + \nabla_z \cdot ( \rho D \frac{\nabla_z \overline{W}}{\overline{W}}) "
              "- Cantera $", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

#%% Soret diffusion
fig, ax = plt.subplots()
vmin = -5E1; vmax = -vmin
im = ax.imshow(T10[:,:], origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$\nabla_x \rho D_T \frac{\nabla_x T}{T} + \nabla_z \rho D_T \frac{\nabla_z T}{T}- Pele$", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')

isp = gas_mix.species_index("H")
fig, ax = plt.subplots()
MW = gas_mix.molecular_weights[isp] / 1E-3
Wbar = gas_mix.mean_molecular_weight / 1E-3
gradT_x = grad_x(fstate.T)
gradT_z = grad_z(fstate.T)
DT = fstate.thermal_diff_coeffs[:,:,isp] * gradT_x / fstate.T
lapT_x = grad_x(DT)
DT = fstate.thermal_diff_coeffs[:,:,isp] * gradT_z / fstate.T
lapT_z = grad_z(DT)

im = ax.imshow(lapT_x + lapT_z, origin="lower",
               vmin=vmin, vmax=vmax,
                cmap="seismic", extent=extent)
ax.set_title(r"$\nabla_x \rho D_T \frac{\nabla_x T}{T} + \nabla_z \rho D_T \frac{\nabla_z T}{T}-Cantera$", fontsize = 18)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')



# %%
